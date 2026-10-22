#!/usr/bin/env python

# Script to filter premapped (STAR) reads before hybrid identification
# A. M. Chakrabarti
# 9th October 2018; 24th December 2020
# October 2026: linear-read filter added; renamed from filter_spliced_reads.py
#
# Each primary STAR record is classified as one of
#   spliced  all of its junctions are annotated (every jM >= 20). Removed, as before.
#   linear   STAR aligned it as a linear read over (nearly) its whole length:
#              no N (intron) operation in the CIGAR,
#              no single deletion longer than --linear-max-del (default 8, what BLAT bridges),
#              aligned query bases (M/I/=/X) >= --linear-min-frac of the read (default 0.9).
#            Such reads are not chimeras. Left in, BLAT can split one carrying a deletion fused
#            with a substitution (a crosslink-type RT error) or a repeat-length change into two
#            abutting hits that pass hybrid identification as a short "loop". Removed.
#   kept     everything else (unmapped, unannotated junctions, partial alignments). Sent to BLAT.
#
# Log (stdout), one count per line:
#   Total reads / Reads kept / Reads discarded (= spliced) / Reads linear
# With --upstream-log (a two-stage run: spliced filter on a genome alignment, then the linear
# filter on a transcriptome alignment of the reads that survived), Total and Reads discarded are
# taken from the upstream log, so the final log still describes the whole sample.

import argparse
import sys
from subprocess import run

import pysam

# pysam CIGAR operation codes
BAM_CMATCH, BAM_CINS, BAM_CDEL, BAM_CREF_SKIP, BAM_CSOFT_CLIP, BAM_CHARD_CLIP, BAM_CPAD, BAM_CEQUAL, BAM_CDIFF = range(9)
ALIGNED_QUERY_OPS = (BAM_CMATCH, BAM_CINS, BAM_CEQUAL, BAM_CDIFF)


def is_annotated_spliced(read):
    """True if every junction in the read is annotated (STAR adds 20 to the motif of an
    annotated junction; -1 means no junction, so unspliced reads are never 'spliced')."""
    if not read.has_tag("jM"):
        return False
    return all(tag >= 20 for tag in read.get_tag("jM"))


def is_linear(read, min_frac, max_del):
    """True if STAR aligned the read linearly over at least min_frac of its length."""
    cigar = read.cigartuples
    if not cigar:
        return False
    if any(op == BAM_CREF_SKIP for op, _ in cigar):
        return False
    deletions = [length for op, length in cigar if op == BAM_CDEL]
    if deletions and max(deletions) > max_del:
        return False
    read_length = read.infer_read_length()  # includes soft and hard clips
    if not read_length:
        return False
    aligned = sum(length for op, length in cigar if op in ALIGNED_QUERY_OPS)
    return aligned / read_length >= min_frac


def has_junction_tag(read):
    """True if the jM tag records at least one junction."""
    return read.has_tag("jM") and any(tag >= 0 for tag in read.get_tag("jM"))


def read_upstream_log(path):
    counts = {}
    with open(path) as f:
        for line in f:
            if ":" in line:
                key, value = line.split(":", 1)
                counts[key.strip()] = int(value.strip())
    return counts


def filter_bam(args):
    stem = args.stem
    bam_in = pysam.AlignmentFile(args.bam, "rb")
    kept_out = pysam.AlignmentFile(f"{stem}.filtered.bam", "wb", template=bam_in)
    spliced_out = pysam.AlignmentFile(f"{stem}.spliced.bam", "wb", template=bam_in) if args.write_discarded else None
    linear_out = pysam.AlignmentFile(f"{stem}.linear.bam", "wb", template=bam_in) if args.write_discarded else None

    total = kept = spliced = linear = 0
    skipped_nonprimary = cigar_tag_disagree = 0

    for read in bam_in.fetch(until_eof=True):

        # STAR is run with --outSAMmultNmax 1, so each read has one primary record. Anything else
        # is not counted, so that read counts stay equal to the trimmed-read counts.
        if read.is_secondary or read.is_supplementary:
            skipped_nonprimary += 1
            continue

        total += 1

        if read.is_unmapped:
            kept += 1
            kept_out.write(read)
            continue

        if read.has_tag("jM"):
            has_n = any(op == BAM_CREF_SKIP for op, _ in (read.cigartuples or []))
            if has_n != has_junction_tag(read):
                cigar_tag_disagree += 1

        if args.spliced_filter and is_annotated_spliced(read):
            spliced += 1
            if spliced_out:
                spliced_out.write(read)
        elif args.linear_filter and is_linear(read, args.linear_min_frac, args.linear_max_del):
            linear += 1
            if linear_out:
                linear_out.write(read)
        else:
            kept += 1
            kept_out.write(read)

    for f in (bam_in, kept_out, spliced_out, linear_out):
        if f:
            f.close()

    if skipped_nonprimary:
        print(f"Warning: skipped {skipped_nonprimary} secondary/supplementary records", file=sys.stderr)
    if cigar_tag_disagree:
        print(f"Warning: {cigar_tag_disagree} records where CIGAR N and the jM tag disagree", file=sys.stderr)

    # A two-stage run reports against the sample's full read count
    if args.upstream_log:
        up = read_upstream_log(args.upstream_log)
        if up["Reads kept"] != total:
            sys.exit(f"Upstream log reports {up['Reads kept']} kept reads but this BAM has {total}")
        report_total = up["Total reads"]
        spliced += up["Reads discarded"]
        linear += up.get("Reads linear", 0)
    else:
        report_total = total

    assert report_total == kept + spliced + linear

    print("Total reads:", report_total)
    print("Reads kept:", kept)
    print("Reads discarded:", spliced)
    print("Reads linear:", linear)

    # Reads kept go on to BLAT as fastq
    run(f"samtools fastq -0 {stem}.filtered.fastq.gz {stem}.filtered.bam", shell=True, check=True)


def main():
    parser = argparse.ArgumentParser(description="Filter premapped (STAR) reads: annotated-spliced and linear reads")
    parser.add_argument("bam", help="STAR BAM (all reads, unmapped included)")
    parser.add_argument("stem", help="Output stem: <stem>.filtered.bam, <stem>.filtered.fastq.gz")
    parser.add_argument("--no-spliced-filter", dest="spliced_filter", action="store_false",
                        help="Do not remove annotated-spliced reads (e.g. transcriptome alignments)")
    parser.add_argument("--no-linear-filter", dest="linear_filter", action="store_false",
                        help="Do not remove linear reads")
    parser.add_argument("--linear-min-frac", type=float, default=0.9,
                        help="Minimum fraction of the read in M/I/=/X for a linear read [0.9]")
    parser.add_argument("--linear-max-del", type=int, default=8,
                        help="Maximum length of any single deletion in a linear read [8]")
    parser.add_argument("--upstream-log", default=None,
                        help="Log of an upstream filtering stage, to report whole-sample totals")
    parser.add_argument("--write-discarded", action="store_true",
                        help="Also write <stem>.spliced.bam and <stem>.linear.bam")
    args = parser.parse_args()
    if not 0 < args.linear_min_frac <= 1:
        parser.error("--linear-min-frac must be in (0, 1]")
    filter_bam(args)


if __name__ == "__main__":
    main()
