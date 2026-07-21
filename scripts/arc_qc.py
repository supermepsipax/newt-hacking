"""Reusable FASTQ QC helpers for Cell Ranger ARC input validation.

These are the functions derived and explained in notebooks/01_arc_input_prep.ipynb, extracted
here so other notebooks and scripts can import them rather than copy them. Notebook 01 remains
the place where the reasoning lives; this is just the library form.
"""

from __future__ import annotations

import gzip
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

ARC_HOME = Path("/home/supermepsipax/Code/cellranger/cellranger-arc-2.2.0")
GEX_WHITELIST = ARC_HOME / "lib/python/cellranger/barcodes/737K-arc-v1.txt.gz"
ATAC_WHITELIST = ARC_HOME / "lib/python/atac/barcodes/737K-arc-v1.txt.gz"

_RC = str.maketrans("ACGTN", "TGCAN")


def revcomp(s: str) -> str:
    return s.translate(_RC)[::-1]


@dataclass(frozen=True)
class ReadSpec:
    length: int
    contents: str
    required: bool = True


GEX_SPEC = {
    "R1": ReadSpec(28, "16 bp cell barcode + 12 bp UMI"),
    "R2": ReadSpec(90, "cDNA insert"),
    "I1": ReadSpec(10, "i7 sample index", required=False),
    "I2": ReadSpec(10, "i5 sample index", required=False),
}
ATAC_SPEC = {
    "R1": ReadSpec(50, "genomic DNA"),
    "R2": ReadSpec(16, "16 bp cell barcode (from i5)"),
    "R3": ReadSpec(49, "genomic DNA (mate of R1)"),
    "I1": ReadSpec(8, "i7 sample index", required=False),
}
SPECS = {"Gene Expression": GEX_SPEC, "Chromatin Accessibility": ATAC_SPEC}

# GEX R1 is 28 bp (barcode + UMI) so the barcode is the first 16; ATAC R2 is the bare barcode.
EXPECTED_INTERP = {"Gene Expression": "first16", "Chromatin Accessibility": "as-is"}

FASTQ_RE = re.compile(
    r"^(?P<sample>.+?)_S(?P<snum>\d+)(?:_L(?P<lane>\d{3}))?"
    r"_(?P<read>[RI][123])_(?P<chunk>\d{3})\.fastq(?P<gz>\.gz)?$"
)


def parse_fastq_name(path) -> dict | None:
    m = FASTQ_RE.match(Path(path).name)
    return m.groupdict() | {"path": str(path)} if m else None


_wl_cache: dict[Path, set[str]] = {}


def load_whitelist(path: Path) -> set[str]:
    if path not in _wl_cache:
        with gzip.open(path, "rt") as fh:
            _wl_cache[path] = {line.strip() for line in fh if line.strip()}
    return _wl_cache[path]


def head_seqs(path: Path, n: int = 5000) -> list[str]:
    """First n sequences. Uses zcat so huge files cost only what we read."""
    proc = subprocess.Popen(f"zcat '{path}' | head -{n * 4}", shell=True,
                            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    seqs = [ln.strip() for i, ln in enumerate(proc.stdout) if i % 4 == 1]
    proc.stdout.close()
    proc.wait()
    return seqs


def barcode_match_rates(path: Path, whitelist: set[str], limit: int = 5000) -> dict[str, float]:
    seqs = head_seqs(path, limit)
    if not seqs:
        return {}
    variants = {
        "as-is": lambda s: s,
        "revcomp": lambda s: revcomp(s),
        "first16": lambda s: s[:16],
        "last16": lambda s: s[-16:],
        "revcomp(first16)": lambda s: revcomp(s[:16]),
        "revcomp(last16)": lambda s: revcomp(s[-16:]),
    }
    return {k: sum(fn(s) in whitelist for s in seqs) / len(seqs) for k, fn in variants.items()}


def preflight(directory, library_type: str, limit: int = 5000) -> bool:
    """Check one library directory against everything Cell Ranger ARC requires."""
    directory = Path(directory)
    spec = SPECS[library_type]
    whitelist = load_whitelist(
        GEX_WHITELIST if library_type == "Gene Expression" else ATAC_WHITELIST)
    bc_read = "R1" if library_type == "Gene Expression" else "R2"
    problems: list[str] = []

    print(f"\n{'=' * 66}\nPREFLIGHT: {directory}  [{library_type}]\n{'=' * 66}")
    files = sorted(directory.glob("*.fastq.gz"))
    if not files:
        print("  FAIL  no .fastq.gz files found")
        return False

    parsed = {f: parse_fastq_name(f) for f in files}
    bad = [f.name for f, p in parsed.items() if p is None]
    if bad:
        problems.append(f"filenames do not match the bcl2fastq pattern: {bad}")
    ok = {f: p for f, p in parsed.items() if p}
    samples = {p["sample"] for p in ok.values()}
    if len(samples) > 1:
        problems.append(f"multiple sample prefixes in one directory: {sorted(samples)}")

    present = {p["read"] for p in ok.values()}
    missing = {r for r, s in spec.items() if s.required} - present
    if missing:
        problems.append(f"missing required reads: {sorted(missing)}")

    for f, p in ok.items():
        exp = spec.get(p["read"])
        if exp is None:
            continue
        seqs = head_seqs(f, 2000)
        if not seqs:
            problems.append(f"{f.name}: no reads")
            continue
        obs = max(set(len(s) for s in seqs), key=[len(s) for s in seqs].count)
        if obs != exp.length:
            note = f"{f.name}: {p['read']} is {obs} bp vs spec {exp.length} bp"
            if obs < exp.length and exp.length <= 28:
                problems.append(note + " -- barcode/UMI reads cannot be short, not repairable")
            else:
                print(f"  note  {note}")

    for f, p in ok.items():
        if p["read"] != bc_read:
            continue
        rates = barcode_match_rates(f, whitelist, limit)
        expected = EXPECTED_INTERP[library_type]
        got = rates[expected]
        best, best_rate = max(rates.items(), key=lambda kv: kv[1])
        print(f"  barcode  {f.name}: {expected!r} (expected) = {got:.1%}; "
              f"best = {best!r} at {best_rate:.1%}")
        if got >= 0.50:
            continue
        if best_rate >= 0.50:
            problems.append(f"{f.name}: barcodes present but only as {best!r}, not the "
                            f"expected {expected!r} -- transform before running the pipeline")
        else:
            problems.append(f"{f.name}: no interpretation matches the {library_type} "
                            f"whitelist (best {best!r} = {best_rate:.1%})")

    if problems:
        print("\n  RESULT: NOT READY")
        for p in problems:
            print(f"    - {p}")
        return False
    print("\n  RESULT: READY")
    return True
