#!/usr/bin/env python
"""Tests for bin/flag_repeat_shifts.py.

Run: python tests/test_flag_repeat_shifts.py   (needs pysam, pandas)
"""

import os
import subprocess
import sys
import tempfile

import pandas as pd
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "bin", "flag_repeat_shifts.py")

# tx1: unique sequence, then a (TG)n repeat (positions 101-140), then unique, then a poly-T run
UNIQUE1 = "ACGATCGTAGCTAGCTTAGCCGATCGATTACGGCATCGATCAGCTAGGCTAACGTTAGCAGTCGATCGTACGATGCATGCTAGTCAGCTAGCATCGACTG"  # 100
REPEAT = "TG" * 20                                                                                                  # 101-140
UNIQUE2 = "CAGTACGATCGGATCAGCATGCAGTACGTAGCATCGATCGTAGCATGCATCAGTCGA"                                                # 141-197
POLYT = "TTTTTTTTTT"                                                                                                # 198-207
UNIQUE3 = "GACTAGCTAGCATCGATCGACTAGCATGCATCGACGTAGCTAGCAGCTAGCATCGATCGAGCTAGCATGCAGT"                               # 208-
TX1 = UNIQUE1 + REPEAT + UNIQUE2 + POLYT + UNIQUE3


def hyb(name, L_start, L_end, R_start, R_end, orientation="genomic", L="tx1", R="tx1"):
    return dict(name=name, orientation=orientation, type="intragenic" if L == R else "intergenic",
                hybrid_selection="single", L_seqnames=L, L_read_start=1, L_read_end=30, L_start=L_start,
                L_end=L_end, L_width=L_end - L_start + 1, L_strand="+", R_seqnames=R, R_read_start=31,
                R_read_end=60, R_start=R_start, R_end=R_end, R_width=R_end - R_start + 1, R_strand="+")


# expected repeat_shift: "TRUE", "FALSE" or "NA" (not tested)
CASES = [
    (hyb("repeat_shift_10", 71, 100, 111, 140), "TRUE"),      # gap 101-110 = TGTGTGTGTG, copies 111-120
    (hyb("repeat_shift_in_repeat", 81, 104, 115, 150), "TRUE"),  # gap 105-114 inside the repeat
    (hyb("homopolymer_1", 170, 199, 201, 230), "TRUE"),        # gap 200 = T, flank T
    (hyb("unique_gap", 1, 40, 46, 80), "FALSE"),               # gap 41-45 in unique sequence
    (hyb("reverse_order", 71, 100, 111, 140, orientation="reverse"), "NA"),
    (hyb("abutting", 71, 100, 101, 140), "NA"),                # gap 0: not tested
    (hyb("long_gap", 1, 30, 60, 90), "NA"),                    # gap 29 > 20
    (hyb("intergenic", 71, 100, 111, 140, R="tx2"), "NA"),
    (hyb("start_clipped", 1, 2, 6, 30), "FALSE"),              # gap 3-5; no upstream flank (clipped), down 6-8 differs
]


def run(hybrids, fasta, out, *extra):
    res = subprocess.run([sys.executable, SCRIPT, hybrids, fasta, out, *extra], capture_output=True, text=True, check=True)
    return dict((k.strip(), int(v)) for k, v in (line.split(":", 1) for line in res.stdout.splitlines()))


def test_all():
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, "tx.fa")
        with open(fa, "w") as f:
            f.write(">tx1\n" + TX1 + "\n>tx2\n" + UNIQUE3 + UNIQUE1 + "\n")
        pysam.faidx(fa)
        hp = os.path.join(d, "h.tsv.gz")
        pd.DataFrame([c for c, _ in CASES]).to_csv(hp, sep="\t", index=False)

        out = os.path.join(d, "flagged.tsv.gz")
        log = run(hp, fa, out)
        got = pd.read_csv(out, sep="\t", keep_default_na=False)
        assert list(got.columns)[-1] == "repeat_shift"
        exp = {c["name"]: e for c, e in CASES}
        assert dict(zip(got["name"], got["repeat_shift"])) == exp, dict(zip(got["name"], got["repeat_shift"]))
        assert log["Hybrids tested"] == sum(e != "NA" for e in exp.values())
        assert log["Repeat shifts flagged"] == sum(e == "TRUE" for e in exp.values())
        assert log["Hybrids out"] == len(CASES)

        out = os.path.join(d, "dropped.tsv.gz")
        log = run(hp, fa, out, "--drop")
        got = pd.read_csv(out, sep="\t", keep_default_na=False)
        assert set(got["name"]) == {n for n, e in exp.items() if e != "TRUE"}
        assert log["Hybrids out"] == len(CASES) - log["Repeat shifts flagged"]

    print("test_flag_repeat_shifts: all passed")


if __name__ == "__main__":
    test_all()
