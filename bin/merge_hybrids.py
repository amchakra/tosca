#!/usr/bin/env python

import gzip
import sys
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', nargs='+', required=True, help='Input .tsv.gz files')
    parser.add_argument('--output', required=True, help='Output .tsv.gz file')
    parser.add_argument('--sample', required=True, help='Sample ID')
    parser.add_argument('--type', required=True, help='Type (e.g., "atlas" to add sample column)')
    return parser.parse_args()

def is_nonempty(file_path):
    return os.path.isfile(file_path) and os.path.getsize(file_path) > 0

def stream_merge(input_files, output_file, sample_id, type_str):
    add_sample = (type_str == "atlas")
    header_written = False

    with gzip.open(output_file, 'wt') as out:
        for file in input_files:
            if not is_nonempty(file):
                sys.stderr.write(f"[merge] Skipping empty or missing file: {file}\n")
                continue
            with gzip.open(file, 'rt') as f:
                for i, line in enumerate(f):
                    line = line.rstrip('\n')
                    if i == 0:
                        if not header_written:
                            if add_sample:
                                out.write(line + '\tsample\n')
                            else:
                                out.write(line + '\n')
                            header_written = True
                        continue  # skip header for remaining files
                    if add_sample:
                        out.write(line + f'\t{sample_id}\n')
                    else:
                        out.write(line + '\n')

def main():
    args = parse_args()
    stream_merge(args.input, args.output, args.sample, args.type)

if __name__ == '__main__':
    main()
