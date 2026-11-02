#!/usr/bin/env python
"""Tests for bin/filter_premap_reads.py.

Run: python tests/test_filter_premap_reads.py   (needs pysam and samtools on PATH)
"""

import array
import gzip
import os
import subprocess
import sys
import tempfile

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "bin", "filter_premap_reads.py")

# name: (cigar or None for unmapped, jM tag, expected class with both filters on)
CASES = {
    "linear_clean":        ("60M",          [-1],      "linear"),
    "linear_1D_mismatch":  ("24M1D26M",     [-1],      "linear"),    # deletion + substitution (MD not needed)
    "linear_2D":           ("42M2D24M",     [-1],      "linear"),
    "linear_softclip_7":   ("7S79M",        [-1],      "linear"),    # 79/86 = 0.92 >= 0.9
    "del_9_kept":          ("41M9D29M",     [-1],      "kept"),      # deletion > 8: a possible 9 nt loop
    "partial_two_arm":     ("30M40S",       [-1],      "kept"),      # one arm aligned: a chimera
    "softclip_15_kept":    ("71M15S",       [-1],      "kept"),      # 71/86 = 0.83 < 0.9
    "annotated_junction":  ("30M500N30M",   [21],      "spliced"),
    "unannotated_N":       ("30M500N30M",   [1],       "kept"),      # novel junction: left for BLAT
    "mixed_junctions":     ("20M300N20M400N20M", [21, 2], "kept"),  # not all annotated
    "unmapped":            (None,           None,      "kept"),
}


def query_length(cigar):
    ops = pysam.AlignedSegment()
    ops.cigarstring = cigar
    return ops.infer_read_length()


def make_bam(path):
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": "chr1", "LN": 100000}]}
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for i, (name, (cigar, jm, _)) in enumerate(CASES.items()):
            seg = pysam.AlignedSegment(bam.header)
            seg.query_name = name
            n = query_length(cigar) if cigar else 60
            seg.query_sequence = "ACGT" * (n // 4) + "ACGT"[: n % 4]
            seg.query_qualities = pysam.qualitystring_to_array("I" * n)
            if cigar is None:
                seg.flag = 4
                seg.reference_id = -1
                seg.reference_start = -1
            else:
                seg.flag = 0
                seg.reference_id = 0
                seg.reference_start = 1000 + 1000 * i
                seg.mapping_quality = 255
                seg.cigarstring = cigar
                seg.set_tag("jM", array.array("b", jm))  # STAR writes jM as B:c
            bam.write(seg)


def run_filter(bam, stem, *extra):
    out = subprocess.run([sys.executable, SCRIPT, bam, stem, *extra], capture_output=True, text=True, check=True)
    counts = {}
    for line in out.stdout.splitlines():
        key, value = line.split(":", 1)
        counts[key.strip()] = int(value)
    return counts


def fastq_names(path):
    with gzip.open(path, "rt") as f:
        return {line[1:].strip() for i, line in enumerate(f) if i % 4 == 0}


def expected(spliced_filter=True, linear_filter=True):
    out = {}
    for name, (_, _, cls) in CASES.items():
        if cls == "spliced" and not spliced_filter:
            cls = "kept"
        if cls == "linear" and not linear_filter:
            cls = "kept"
        out[name] = cls
    return out


def check(counts, fastq, exp, total=len(CASES)):
    kept = {n for n, c in exp.items() if c == "kept"}
    assert counts["Total reads"] == total, counts
    assert counts["Reads kept"] == len(kept), counts
    assert counts["Reads discarded"] == sum(c == "spliced" for c in exp.values()), counts
    assert counts["Reads linear"] == sum(c == "linear" for c in exp.values()), counts
    assert counts["Reads kept"] + counts["Reads discarded"] + counts["Reads linear"] == counts["Total reads"]
    assert fastq_names(fastq) == kept, (fastq_names(fastq), kept)


def test_all():
    with tempfile.TemporaryDirectory() as d:
        bam = os.path.join(d, "in.bam")
        make_bam(bam)

        # Both filters (default)
        stem = os.path.join(d, "both")
        counts = run_filter(bam, stem, "--write-discarded")
        check(counts, stem + ".filtered.fastq.gz", expected())
        linear_names = {r.query_name for r in pysam.AlignmentFile(stem + ".linear.bam", check_sq=False)}
        assert linear_names == {n for n, c in expected().items() if c == "linear"}
        spliced_names = {r.query_name for r in pysam.AlignmentFile(stem + ".spliced.bam", check_sq=False)}
        assert spliced_names == {"annotated_junction"}

        # Previous behaviour: spliced filter only
        stem = os.path.join(d, "spliced_only")
        check(run_filter(bam, stem, "--no-linear-filter"), stem + ".filtered.fastq.gz", expected(linear_filter=False))

        # Transcriptome path: linear filter only
        stem = os.path.join(d, "linear_only")
        check(run_filter(bam, stem, "--no-spliced-filter"), stem + ".filtered.fastq.gz", expected(spliced_filter=False))

        # Thresholds: --linear-max-del 9 makes the 9 nt deletion linear; --linear-min-frac 0.8 the 15 nt clip
        stem = os.path.join(d, "loose")
        exp = expected()
        exp["del_9_kept"] = "linear"
        exp["softclip_15_kept"] = "linear"
        check(run_filter(bam, stem, "--linear-max-del", "9", "--linear-min-frac", "0.8"), stem + ".filtered.fastq.gz", exp)

        # Two-stage run: stage 1 spliced filter, stage 2 linear filter on what survived, with the
        # stage-1 log giving the whole-sample totals
        stem1 = os.path.join(d, "stage1")
        c1 = run_filter(bam, stem1, "--no-linear-filter")
        log1 = stem1 + ".log"
        with open(log1, "w") as f:
            for k, v in c1.items():
                f.write(f"{k}: {v}\n")
        stage1_bam = stem1 + ".filtered.bam"
        stem2 = os.path.join(d, "stage2")
        c2 = run_filter(stage1_bam, stem2, "--no-spliced-filter", "--upstream-log", log1)
        check(c2, stem2 + ".filtered.fastq.gz", expected())

    print("test_filter_premap_reads: all passed")


if __name__ == "__main__":
    test_all()
