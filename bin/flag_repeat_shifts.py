#!/usr/bin/env python

# Script to flag repeat-shift hybrids
# A. M. Chakrabarti
# October 2026
#
# A linear read through a tandem repeat that carries fewer repeat copies than the reference
# (repeat-length variant, polymerase slippage) can come out of BLAT as two abutting blocks: the
# upstream block stops where the repeat starts and the downstream block carries the read's repeat,
# placed on the 3' end of the reference repeat. Hybrid identification then reports a
# genomic-order, same-transcript hybrid whose "loop" is the missing repeat copies.
#
# Tell-tale: the k skipped transcript bases (k = transcript gap) are an exact copy of the k bases
# immediately upstream or downstream of them. Tested only for intramolecular, genomic-order
# hybrids with 1 <= gap <= --max-gap (default 20); repeat_shift is TRUE/FALSE for those and NA
# otherwise. A real hairpin loop matches its flank only by chance (reverse-order hybrids, which
# cannot come from a linear read, show ~0% at gaps > 3 nt).
#
# Adds a repeat_shift column (last). With --drop, flagged hybrids are removed instead.

import argparse
import sys

import numpy as np
import pandas as pd
import pysam


def flag_repeat_shifts(hybrids, fasta, max_gap):
    flags = pd.Series(pd.NA, index=hybrids.index, dtype="object")
    if hybrids.empty:
        return flags

    upstream_end = np.minimum(hybrids["L_end"], hybrids["R_end"])           # 1-based end of the upstream arm
    downstream_start = np.maximum(hybrids["L_start"], hybrids["R_start"])
    tx_gap = downstream_start - upstream_end - 1

    tested = ((hybrids["L_seqnames"] == hybrids["R_seqnames"]) &
              (hybrids["orientation"] == "genomic") &
              (tx_gap >= 1) & (tx_gap <= max_gap))

    lengths = dict(zip(fasta.references, fasta.lengths))
    for i in hybrids.index[tested]:
        seqname = hybrids.at[i, "L_seqnames"]
        g0 = int(upstream_end[i])        # 0-based start of the gap
        k = int(tx_gap[i])
        n = lengths[seqname]
        gap = fasta.fetch(seqname, g0, g0 + k).upper()
        up = fasta.fetch(seqname, g0 - k, g0).upper() if g0 - k >= 0 else None
        down = fasta.fetch(seqname, g0 + k, g0 + 2 * k).upper() if g0 + 2 * k <= n else None
        flags.at[i] = (up is not None and gap == up) or (down is not None and gap == down)

    return flags


def main():
    parser = argparse.ArgumentParser(description="Flag genomic-order hybrids whose gap is a tandem-repeat copy")
    parser.add_argument("hybrids", help="Hybrids TSV(.gz) with transcript coordinates")
    parser.add_argument("fasta", help="Transcript FASTA (indexed; .fai alongside)")
    parser.add_argument("output", help="Output TSV(.gz)")
    parser.add_argument("--max-gap", type=int, default=20, help="Largest transcript gap tested [20]")
    parser.add_argument("--drop", action="store_true", help="Remove flagged hybrids instead of only flagging")
    args = parser.parse_args()

    hybrids = pd.read_csv(args.hybrids, sep="\t")
    fasta = pysam.FastaFile(args.fasta)
    flags = flag_repeat_shifts(hybrids, fasta, args.max_gap)

    tested = int(flags.notna().sum())
    flagged = int((flags == True).sum())  # noqa: E712 (object dtype with NA)

    hybrids["repeat_shift"] = flags.map({True: "TRUE", False: "FALSE"}).fillna("NA")
    if args.drop:
        hybrids = hybrids[hybrids["repeat_shift"] != "TRUE"]

    hybrids.to_csv(args.output, sep="\t", index=False)

    print(f"Hybrids in: {len(flags)}")
    print(f"Hybrids tested: {tested}")
    print(f"Repeat shifts flagged: {flagged}")
    print(f"Repeat shifts dropped: {flagged if args.drop else 0}")
    print(f"Hybrids out: {len(hybrids)}")


if __name__ == "__main__":
    sys.exit(main())
