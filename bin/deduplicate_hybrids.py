#!/usr/bin/env python

# Script to deduplcate using UMI tools approach adapted to hybrids
# A. M. Chakrabarti
# 25th December 2020

import sys
import re
import pandas as pd
import time
from umi_tools import UMIClusterer

def deduplicate(grouped_hybrids):
    """Deduplicate group of hybrids with same mapping coordinates"""

    # Get UMI dictionary with counts and convert to bytes for UMI tools
    umi_d = grouped_hybrids['umi'].value_counts().to_dict()
    umi_d = {k.encode('utf-8'):v for k,v in umi_d.items()}

    # Cluster UMIs
    dedup = clusterer(umi_d, threshold=1)

    # Select representative read per UMI
    dedup_l = []
    for umi in dedup:
        umi_g = [x.decode('utf-8') for x in umi]
        read = grouped_hybrids[grouped_hybrids['umi'].isin(umi_g)].iloc[0]
        dedup_l.append(read)

    dedup_pd = pd.DataFrame(dedup_l) 
    return dedup_pd

# ==========
# Run
# ==========

if len(sys.argv) == 5:

    f_in = sys.argv[1]
    f_out = sys.argv[2]
    umi_separator = sys.argv[3]
    dedup_method = sys.argv[4]

    tic = time.perf_counter()

    # Read in hybrids
    hybrids = pd.read_csv(f_in, sep='\t')
    hybrids = hybrids[hybrids['hybrid_selection'].isin(['single', 'multi_overlap'])]

    if dedup_method != 'none':

        hybrids["umi"] = hybrids['name'].str.replace(".*" + umi_separator, "", regex = True)

        # hybrids_grp = hybrids.groupby(['L_seqnames', 'L_start', 'L_end', 'R_seqnames', 'R_start', 'R_end'])
        # Removed end as requirement as may have been single end sequenced to different lengths
        # Or 3' quality trimmed slightly differently
        hybrids_grp = hybrids.groupby(['L_seqnames', 'L_start', 'R_seqnames', 'R_start'])

        # Deduplicate
        clusterer = UMIClusterer(cluster_method=dedup_method)
        unique_hybrids = [deduplicate(group) for name, group in hybrids_grp]
        unique_hybrids = pd.concat(unique_hybrids)

    else:

        unique_hybrids = hybrids

    toc = time.perf_counter()

    # Write out
    unique_hybrids.to_csv(f_out, sep = '\t', index = False)

    print(f"Deduplicated in {toc - tic:0.2f} seconds")
    print(f"Hybrids in: {hybrids.shape[0]}")
    print(f"Hybrids out: {unique_hybrids.shape[0]}")
    print(f"PCR duplication ratio: {hybrids.shape[0]/unique_hybrids.shape[0]:0.2f}")

else:

    print("deduplicate_hybrids.py <input_hybrids> <output_hybrids> <umi_separator> <deduplication_method>")