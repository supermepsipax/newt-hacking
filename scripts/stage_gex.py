#!/usr/bin/env python3
"""Stage the MGI GEX (RNA) FASTQs into the layout Cell Ranger ARC expects.

Unlike the ATAC side, the RNA reads need **no transformation at all** -- notebook 02 §9
measured them at 28 bp / 90 bp with 88.7% of R1 barcodes matching the GEX whitelist in the
forward orientation, and Cell Ranger ignores MGI's trailing "/1" and "/2" (see
`split_atac.header_prefix`). So the only thing standing between the delivery and the pipeline
is the *filename*:

    SI-TT-G9_L01_R1.fastq.gz   ->   5wks_f_GEX_S1_L001_R1_001.fastq.gz

which means symlinks are enough. Rewriting the bytes would cost a full decompress/recompress
pass over ~100 GB to produce identical data.

Two things this does that a shell one-liner tends to get wrong:

* **Lane comes from the filename, not the directory.** MGI writes `_L01_` where bcl2fastq
  wants `_L001_`. Collapsing four lanes onto one name silently discards three quarters of the
  library, and nothing downstream would flag it.
* **It validates before it commits.** Read lengths and whitelist match rate are checked on a
  sample of reads, so a swapped index or a mislabelled file surfaces here rather than three
  hours into `cellranger-arc count`.

Usage:
    python scripts/stage_gex.py --in-dir  /path/to/MGI_Outs/RNA \
                                --out-dir /path/to/arc_ready/gex \
                                --map SI-TT-G9=5wks_f --map SI-TT-B9=9wks_f

Add --dry-run to see the plan without touching the filesystem, or --copy to materialise real
files instead of symlinks (only needed if the destination will be read somewhere the source
path is not mounted).
"""

from __future__ import annotations

import argparse
import gzip
import re
import subprocess
import sys
from pathlib import Path

ARC_HOME = Path(__file__).resolve().parent.parent / "cellranger-arc-2.2.0"
GEX_WHITELIST = ARC_HOME / "lib/python/cellranger/barcodes/737K-arc-v1.txt.gz"

# SI-TT-G9_L01_R1.fastq.gz -> index, lane, read
NAME_RE = re.compile(r"^(?P<idx>.+?)_L0*(?P<lane>\d{1,3})_(?P<read>R[12])\.fastq\.gz$")

EXPECTED_LEN = {"R1": 28, "R2": 90}      # 16 bp barcode + 12 bp UMI; cDNA insert
MIN_BC_MATCH = 0.50


def load_whitelist() -> set[str]:
    with gzip.open(GEX_WHITELIST, "rt") as fh:
        return {line.strip() for line in fh if line.strip()}


def head_seqs(path: Path, n: int = 2000) -> list[str]:
    """First n sequences. Streams, so a 40 GB file costs the same as a 40 MB one."""
    proc = subprocess.Popen(f"zcat '{path}' | head -{n * 4}", shell=True,
                            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    seqs = [ln.strip() for i, ln in enumerate(proc.stdout) if i % 4 == 1]
    proc.stdout.close()
    proc.wait()
    return seqs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-dir", required=True, type=Path,
                    help="RNA delivery root; searched recursively for *.fastq.gz")
    ap.add_argument("--out-dir", required=True, type=Path,
                    help="staging root; one subdirectory is created per sample")
    ap.add_argument("--map", required=True, action="append", metavar="INDEX=SAMPLE",
                    help="plate index to biological sample, e.g. SI-TT-G9=5wks_f (repeatable)")
    ap.add_argument("--suffix", default="_GEX",
                    help="appended to the sample name to form the FASTQ prefix")
    ap.add_argument("--copy", action="store_true",
                    help="copy the files instead of symlinking them")
    ap.add_argument("--dry-run", action="store_true", help="print the plan and stop")
    ap.add_argument("--check-reads", type=int, default=2000)
    args = ap.parse_args()

    try:
        mapping = dict(m.split("=", 1) for m in args.map)
    except ValueError:
        print("ERROR: --map takes INDEX=SAMPLE, e.g. --map SI-TT-G9=5wks_f", file=sys.stderr)
        return 1

    files = sorted(args.in_dir.rglob("*.fastq.gz"))
    if not files:
        print(f"ERROR: no .fastq.gz under {args.in_dir}", file=sys.stderr)
        return 1

    # ---- build the plan -----------------------------------------------------
    plan: list[tuple[Path, Path, str]] = []      # (src, dst, read)
    unmapped: set[str] = set()
    for f in files:
        m = NAME_RE.match(f.name)
        if not m:
            print(f"  skip  {f.name}  (does not look like INDEX_Lnn_R#.fastq.gz)")
            continue
        idx, lane, read = m["idx"], int(m["lane"]), m["read"]
        sample = mapping.get(idx)
        if sample is None:
            unmapped.add(idx)
            continue
        prefix = f"{sample}{args.suffix}"
        dst = args.out_dir / sample / f"{prefix}_S1_L{lane:03d}_{read}_001.fastq.gz"
        plan.append((f, dst, read))

    if unmapped:
        print(f"ERROR: no --map entry for index(es): {sorted(unmapped)}", file=sys.stderr)
        return 1
    if not plan:
        print("ERROR: nothing to stage", file=sys.stderr)
        return 1

    collisions = {d for _s, d, _r in plan if [d for _s2, d2, _r2 in plan if d2 == d].count(d) > 1}
    if collisions:
        print(f"ERROR: two inputs map to the same output name: {sorted(collisions)}",
              file=sys.stderr)
        return 1

    print(f"{len(plan)} file(s) to stage "
          f"({'copy' if args.copy else 'symlink'}{', DRY RUN' if args.dry_run else ''}):\n")
    for src, dst, _read in plan:
        print(f"  {src.name:<28} -> {dst.relative_to(args.out_dir)}")

    # ---- validate before committing ----------------------------------------
    print("\nvalidating (lengths, and barcodes against the GEX whitelist) ...")
    whitelist = load_whitelist()
    problems: list[str] = []
    for src, dst, read in plan:
        seqs = head_seqs(src, args.check_reads)
        if not seqs:
            problems.append(f"{src.name}: no reads")
            continue
        obs = max(set(len(s) for s in seqs), key=[len(s) for s in seqs].count)
        exp = EXPECTED_LEN[read]
        if obs != exp:
            problems.append(f"{src.name}: {read} is {obs} bp, expected {exp} bp")
        if read == "R1":
            rate = sum(s[:16] in whitelist for s in seqs) / len(seqs)
            flag = "ok" if rate >= MIN_BC_MATCH else "FAIL"
            print(f"  {src.name:<28} {read} {obs:>3} bp   barcode {rate:6.1%}  {flag}")
            if rate < MIN_BC_MATCH:
                problems.append(f"{src.name}: only {rate:.1%} of first-16 barcodes are in the "
                                "GEX whitelist -- wrong file, wrong whitelist, or wrong "
                                "orientation")
        else:
            print(f"  {src.name:<28} {read} {obs:>3} bp")

    if problems:
        print("\nNOT STAGED -- fix these first:", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        return 1

    if args.dry_run:
        print("\ndry run: validation passed, nothing written")
        return 0

    # ---- commit -------------------------------------------------------------
    for src, dst, _read in plan:
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        if args.copy:
            dst.write_bytes(src.read_bytes())
        else:
            dst.symlink_to(src.resolve())

    samples = sorted({d.parent for _s, d, _r in plan})
    print(f"\nstaged {len(plan)} file(s) into {len(samples)} sample director(ies):")
    for d in samples:
        print(f"  {d}  ({len(list(d.glob('*.fastq.gz')))} files)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
