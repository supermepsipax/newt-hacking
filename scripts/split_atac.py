#!/usr/bin/env python3
"""Split MGI 2-file ATAC output into the 3 files Cell Ranger ARC expects.

The MGI/mgikit delivery packs three pieces of information into read 2:

    R2 (73 bp) = [0:49] genomic mate | [49:57] CAGACGCG spacer | [57:73] barcode (revcomp)

Cell Ranger ARC instead wants:

    R1 (50 bp) genomic
    R2 (16 bp) cell barcode, forward orientation
    R3 (49 bp) genomic mate

So this script reads the two input files in lockstep and writes three, renaming them to the
bcl2fastq convention on the way out. It validates barcodes against the shipped whitelist as it
streams, so a mistake surfaces in the first few seconds rather than after several hours.

Usage:
    python scripts/split_atac.py --r1 IN_R1.fastq.gz --r2 IN_R2.fastq.gz \
        --out-dir /path/to/out --sample SI-NA-E10

See notebooks/02_atac_split.ipynb for the derivation of the layout and the correctness tests.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------- layout constants
# Derived empirically in notebook 02 -- do not change without re-running that analysis.
GENOMIC_END = 49       # R2[0:49]   -> new R3
SPACER = b"CAGACGCG"   # R2[49:57]  -> constant, discarded
SPACER_START, SPACER_END = 49, 57
BC_START, BC_END = 57, 73   # R2[57:73] -> revcomp -> new R2
EXPECTED_R1_LEN = 50
EXPECTED_R2_LEN = 73

COMPLEMENT = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")

ARC_HOME = Path("/home/supermepsipax/Code/cellranger/cellranger-arc-2.2.0")
ATAC_WHITELIST = ARC_HOME / "lib/python/atac/barcodes/737K-arc-v1.txt.gz"


def revcomp(seq: bytes) -> bytes:
    return seq.translate(COMPLEMENT)[::-1]


def strip_suffix(header: bytes) -> bytes:
    """Drop the trailing /1 or /2 so all three files share an identical read name.

    Cell Ranger compares `name.split()[0]` across reads (tenkit/fasta.py) and raises
    FastqParseError if they differ, so '/1' vs '/2' vs '/3' would be rejected.
    """
    if len(header) > 2 and header[-2] == 0x2F:  # ord('/')
        return header[:-2]
    return header


def load_whitelist() -> set[bytes]:
    with gzip.open(ATAC_WHITELIST, "rb") as fh:
        return {line.strip() for line in fh if line.strip()}


def open_writer(path: Path, level: int, threads: int):
    """Return (process, file_handle) writing gzip to `path`. Uses pigz when available."""
    fh = open(path, "wb")
    if shutil.which("pigz"):
        cmd = ["pigz", f"-{level}", "-p", str(threads)]
    else:
        cmd = ["gzip", f"-{level}"]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=fh)
    return proc, fh


def open_reader(path: Path):
    """Stream-decompress `path`. A separate zcat process keeps decompression off our core."""
    cmd = ["unpigz", "-c"] if shutil.which("unpigz") else ["zcat"]
    proc = subprocess.Popen([*cmd, str(path)], stdout=subprocess.PIPE, bufsize=1 << 20)
    return proc


def fastq_records(stream):
    """Yield (header, seq, qual) as bytes with trailing newlines stripped."""
    it = iter(stream)
    for header in it:
        seq = next(it)
        next(it)  # '+'
        qual = next(it)
        yield header.rstrip(b"\n"), seq.rstrip(b"\n"), qual.rstrip(b"\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--r1", required=True, type=Path, help="input R1 (50 bp genomic)")
    ap.add_argument("--r2", required=True, type=Path, help="input R2 (73 bp composite)")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--sample", required=True, help="sample prefix for output filenames")
    ap.add_argument("--lane", default=1, type=int, help="lane number for the L00N field")
    ap.add_argument("--limit", type=int, default=None, help="stop after N reads (for testing)")
    ap.add_argument("--gzip-level", type=int, default=3,
                    help="1=fastest, 6=default. 3 is ~4x faster than 6 for ~10%% more size")
    ap.add_argument("--threads", type=int, default=2, help="threads per pigz process")
    ap.add_argument("--check-reads", type=int, default=100_000,
                    help="validate barcodes against the whitelist over the first N reads")
    ap.add_argument("--min-match", type=float, default=0.50,
                    help="abort if whitelist match rate over the check window is below this")
    ap.add_argument("--progress-every", type=int, default=5_000_000)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{args.sample}_S1_L{args.lane:03d}"
    out_paths = {r: args.out_dir / f"{stem}_{r}_001.fastq.gz" for r in ("R1", "R2", "R3")}

    for r, p in out_paths.items():
        if p.exists():
            print(f"ERROR: {p} already exists -- refusing to overwrite.", file=sys.stderr)
            return 1

    print(f"[{time.strftime('%H:%M:%S')}] loading ATAC whitelist ...", file=sys.stderr, flush=True)
    whitelist = load_whitelist()

    print(f"[{time.strftime('%H:%M:%S')}] in  R1 {args.r1}", file=sys.stderr)
    print(f"[{time.strftime('%H:%M:%S')}] in  R2 {args.r2}", file=sys.stderr)
    for r, p in out_paths.items():
        print(f"[{time.strftime('%H:%M:%S')}] out {r} {p}", file=sys.stderr)
    print(f"[{time.strftime('%H:%M:%S')}] compressor: "
          f"{'pigz' if shutil.which('pigz') else 'gzip'} -{args.gzip_level}",
          file=sys.stderr, flush=True)

    rp1, rp2 = open_reader(args.r1), open_reader(args.r2)
    writers = {r: open_writer(p, args.gzip_level, args.threads) for r, p in out_paths.items()}
    w1, w2, w3 = (writers[r][0].stdin for r in ("R1", "R2", "R3"))

    n = 0
    bc_hits = 0
    bc_checked = 0
    spacer_hits = 0
    bad_len = 0
    t0 = time.time()
    aborted = None

    try:
        for (h1, s1, q1), (h2, s2, q2) in zip(fastq_records(rp1.stdout),
                                              fastq_records(rp2.stdout)):
            name = strip_suffix(h1)
            if name != strip_suffix(h2):
                aborted = (f"read {n}: R1/R2 names diverge -- "
                           f"{h1.decode(errors='replace')} vs {h2.decode(errors='replace')}")
                break
            if len(s1) != EXPECTED_R1_LEN or len(s2) != EXPECTED_R2_LEN:
                bad_len += 1
                if bad_len <= 3:
                    print(f"  warn read {n}: lengths R1={len(s1)} R2={len(s2)} "
                          f"(expected {EXPECTED_R1_LEN}/{EXPECTED_R2_LEN})",
                          file=sys.stderr, flush=True)

            barcode = revcomp(s2[BC_START:BC_END])
            bc_qual = q2[BC_START:BC_END][::-1]   # reversed to stay aligned with the revcomp

            w1.write(b"%s\n%s\n+\n%s\n" % (name, s1, q1))
            w2.write(b"%s\n%s\n+\n%s\n" % (name, barcode, bc_qual))
            w3.write(b"%s\n%s\n+\n%s\n" % (name, s2[:GENOMIC_END], q2[:GENOMIC_END]))

            if bc_checked < args.check_reads:
                bc_checked += 1
                bc_hits += barcode in whitelist
                spacer_hits += s2[SPACER_START:SPACER_END] == SPACER
                if bc_checked == args.check_reads:
                    rate = bc_hits / bc_checked
                    sp = spacer_hits / bc_checked
                    print(f"[{time.strftime('%H:%M:%S')}] validation over first "
                          f"{bc_checked:,} reads: barcode whitelist {rate:.1%}, "
                          f"spacer {sp:.1%}", file=sys.stderr, flush=True)
                    if rate < args.min_match:
                        aborted = (f"barcode match rate {rate:.1%} is below "
                                   f"--min-match {args.min_match:.0%}; the layout constants "
                                   "are probably wrong for this file")
                        break

            n += 1
            if args.limit and n >= args.limit:
                break
            if n % args.progress_every == 0:
                el = time.time() - t0
                print(f"[{time.strftime('%H:%M:%S')}] {n:,} reads  "
                      f"{n/el/1000:.0f}k reads/s  elapsed {el/60:.1f} min",
                      file=sys.stderr, flush=True)
    finally:
        for proc, fh in writers.values():
            if proc.stdin:
                proc.stdin.close()
            proc.wait()
            fh.close()
        for rp in (rp1, rp2):
            if rp.stdout:
                rp.stdout.close()
            rp.terminate()

    el = time.time() - t0
    if aborted:
        print(f"\nABORTED: {aborted}", file=sys.stderr)
        print("Partial outputs left in place for inspection; delete before re-running.",
              file=sys.stderr)
        return 1

    print(f"\n[{time.strftime('%H:%M:%S')}] done: {n:,} reads in {el/60:.1f} min "
          f"({n/el/1000:.0f}k reads/s)", file=sys.stderr)
    if bc_checked:
        print(f"  barcode whitelist match (first {bc_checked:,}): "
              f"{bc_hits/bc_checked:.1%}", file=sys.stderr)
    if bad_len:
        print(f"  WARNING: {bad_len:,} reads had unexpected lengths", file=sys.stderr)
    for r, p in out_paths.items():
        print(f"  {r}: {p}  ({p.stat().st_size/1e9:.1f} GB)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
