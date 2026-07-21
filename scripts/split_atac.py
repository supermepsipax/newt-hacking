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

Read names are written through verbatim, MGI's trailing "/1" and "/2" included: Cell Ranger
compares only the part of the name before the first space *or slash*, so the suffixes are
already invisible to it. See `header_prefix` below.

Usage:
    python scripts/split_atac.py --r1 IN_R1.fastq.gz --r2 IN_R2.fastq.gz \
        --out-dir /path/to/out --sample SI-NA-E10

See notebooks/02_atac_split.ipynb for the derivation of the layout and the correctness tests.
"""

from __future__ import annotations

import argparse
import gzip
import re
import shutil
import subprocess
import sys
import time
from itertools import zip_longest
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

ARC_HOME = Path(__file__).resolve().parent.parent / "cellranger-arc-2.2.0"
ATAC_WHITELIST = ARC_HOME / "lib/python/atac/barcodes/737K-arc-v1.txt.gz"

LANE_RE = re.compile(r"_L0*(\d{1,3})[_.]")


def revcomp(seq: bytes) -> bytes:
    return seq.translate(COMPLEMENT)[::-1]


def header_prefix(header: bytes) -> bytes:
    """The part of a read name Cell Ranger uses to match records across files.

    Mirrors 10x's `fastq_set::read_pair_iter`, which does exactly:

        header.split(|x| *x == b' ' || *x == b'/').next()

    That splits on '/' as well as space, so MGI's '@NAME/1' and '@NAME/2' already
    compare equal and need no rewriting -- support for this landed in fastq_set in
    Aug 2020 ("tolerate a slash as a delimiter between prefix and suffix of FASTQ
    header"), well before cellranger-arc 2.0.2. We use it only to *check* that the
    two inputs are in lockstep; the headers themselves are written through verbatim.
    """
    for i, ch in enumerate(header):
        if ch in (0x20, 0x2F):  # ord(' '), ord('/')
            return header[:i]
    return header


def lane_from_filename(path: Path, fallback: int = 1) -> int:
    """Recover the lane number from an MGI filename, e.g. SI-NA-E10_L01_R1 -> 1.

    MGI writes _L01_ where bcl2fastq wants _L001_, and a run spanning several lanes
    must keep them distinct or Cell Ranger sees one library as many duplicates of the
    same lane. Returns `fallback` if the name carries no lane field.
    """
    m = LANE_RE.search(path.name)
    return int(m.group(1)) if m else fallback


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


class TruncatedFastq(Exception):
    """A FASTQ stream ended part-way through a 4-line record."""


def fastq_records(stream, which: str):
    """Yield (header, seq, qual) as bytes with trailing newlines stripped.

    A stream that dies mid-record (drive disconnect, corrupt gzip member) leaves a partial
    entry. Letting the bare `next()` raise StopIteration here would surface as PEP 479's
    opaque "generator raised StopIteration" RuntimeError, so we name the failure instead.
    """
    it = iter(stream)
    for header in it:
        try:
            seq = next(it)
            next(it)  # '+'
            qual = next(it)
        except StopIteration:
            raise TruncatedFastq(
                f"input {which} ended part-way through a FASTQ record -- the stream was cut "
                "off mid-entry"
            ) from None
        yield header.rstrip(b"\n"), seq.rstrip(b"\n"), qual.rstrip(b"\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--r1", required=True, type=Path, help="input R1 (50 bp genomic)")
    ap.add_argument("--r2", required=True, type=Path, help="input R2 (73 bp composite)")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--sample", required=True, help="sample prefix for output filenames")
    ap.add_argument("--lane", default=None, type=int,
                    help="lane number for the L00N field; default is parsed from the input "
                         "filename (_L01_ -> 1), falling back to 1")
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
    lane = args.lane if args.lane is not None else lane_from_filename(args.r1)
    if args.lane is None:
        print(f"[{time.strftime('%H:%M:%S')}] lane {lane} (parsed from {args.r1.name})",
              file=sys.stderr)
    stem = f"{args.sample}_S1_L{lane:03d}"
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
    stopped_early = False   # True only when we chose to stop (--limit or an abort)

    try:
        # zip_longest, not zip: plain zip() ends silently when either input runs out, which
        # turns a mid-run decompressor crash (e.g. the drive dropping off the USB bus) into a
        # clean-looking "done". The None fill lets us tell truncation from normal completion.
        for rec1, rec2 in zip_longest(fastq_records(rp1.stdout, "R1"),
                                      fastq_records(rp2.stdout, "R2"), fillvalue=None):
            if rec1 is None or rec2 is None:
                short = "R1" if rec1 is None else "R2"
                aborted = (f"input {short} ran out after {n:,} reads while the other still had "
                           "records -- the inputs are truncated or mismatched")
                stopped_early = True
                break
            (h1, s1, q1), (h2, s2, q2) = rec1, rec2
            if header_prefix(h1) != header_prefix(h2):
                aborted = (f"read {n}: R1/R2 names diverge -- "
                           f"{h1.decode(errors='replace')} vs {h2.decode(errors='replace')}")
                stopped_early = True
                break
            if len(s1) != EXPECTED_R1_LEN or len(s2) != EXPECTED_R2_LEN:
                bad_len += 1
                if bad_len <= 3:
                    print(f"  warn read {n}: lengths R1={len(s1)} R2={len(s2)} "
                          f"(expected {EXPECTED_R1_LEN}/{EXPECTED_R2_LEN})",
                          file=sys.stderr, flush=True)

            barcode = revcomp(s2[BC_START:BC_END])
            bc_qual = q2[BC_START:BC_END][::-1]   # reversed to stay aligned with the revcomp

            # Headers pass through untouched: R1 keeps its own, the two reads carved out
            # of input R2 keep R2's. This mirrors real 10x output, where each file carries
            # its own read-number suffix (`1:N:0:0` vs `4:N:0:0`) rather than a shared name.
            w1.write(b"%s\n%s\n+\n%s\n" % (h1, s1, q1))
            w2.write(b"%s\n%s\n+\n%s\n" % (h2, barcode, bc_qual))
            w3.write(b"%s\n%s\n+\n%s\n" % (h2, s2[:GENOMIC_END], q2[:GENOMIC_END]))

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
                        stopped_early = True
                        break

            n += 1
            if args.limit and n >= args.limit:
                stopped_early = True
                break
            if n % args.progress_every == 0:
                el = time.time() - t0
                print(f"[{time.strftime('%H:%M:%S')}] {n:,} reads  "
                      f"{n/el/1000:.0f}k reads/s  elapsed {el/60:.1f} min",
                      file=sys.stderr, flush=True)
    except TruncatedFastq as exc:
        aborted = str(exc)
        stopped_early = True
    finally:
        writer_rc = {}
        for r, (proc, fh) in writers.items():
            if proc.stdin:
                proc.stdin.close()
            writer_rc[r] = proc.wait()
            fh.close()
        reader_rc = {}
        for name, rp in (("R1", rp1), ("R2", rp2)):
            if rp.stdout:
                rp.stdout.close()
            if stopped_early:
                rp.terminate()      # we quit on purpose; SIGPIPE/kill is expected
                rp.wait()
            else:
                # We consumed both inputs to EOF, so each decompressor should have exited 0.
                # Anything else means it died mid-stream and our output is truncated.
                reader_rc[name] = rp.wait()

    # A decompressor that dies (I/O error, drive disconnect, corrupt member) closes its stdout,
    # which looks exactly like end-of-input. The exit code is the only way to tell them apart --
    # without this check a half-converted file reports success. This is not hypothetical: it
    # happened when the host suspended and took the USB drive down with it.
    for name, rc in reader_rc.items():
        if rc != 0 and not aborted:
            aborted = (f"decompressor for {name} exited {rc} -- the input was truncated or "
                       "unreadable mid-stream (check `dmesg` for I/O errors, and whether the "
                       "drive dropped off the bus). THE OUTPUT IS INCOMPLETE.")
    for r, rc in writer_rc.items():
        if rc != 0 and not aborted:
            aborted = (f"compressor for output {r} exited {rc} -- the output is incomplete "
                       "(out of disk space, or the destination went away).")

    el = time.time() - t0
    if aborted:
        print(f"\nABORTED after {n:,} reads in {el/60:.1f} min: {aborted}", file=sys.stderr)
        print("The partial outputs below are NOT usable -- delete them before re-running:",
              file=sys.stderr)
        for r, p in out_paths.items():
            if p.exists():
                print(f"  {p}  ({p.stat().st_size/1e9:.1f} GB)", file=sys.stderr)
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
