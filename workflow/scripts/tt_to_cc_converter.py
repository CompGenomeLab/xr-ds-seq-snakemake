#!/usr/bin/env python3
"""
Convert every occurrence of the dinucleotide TT to CC in FASTA sequences.
For each read, generates one variant per TT position, plus the original sequence.

Usage:
    tt_to_cc_converter.py <input.fasta[.gz]> <output.fasta.gz>
"""

import argparse
import gzip
from pathlib import Path
import sys


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t")
    return open(path, mode, encoding="utf-8")


def parse_args():
    parser = argparse.ArgumentParser(description="Convert TT→CC in FASTA reads.")
    parser.add_argument("input_fasta", type=Path, help="Input FASTA (optionally gzipped).")
    parser.add_argument("output_fasta", type=Path, help="Output FASTA (gzipped recommended).")
    return parser.parse_args()


def main():
    args = parse_args()

    infile = args.input_fasta
    outfile = args.output_fasta

    with open_text(infile, "r") as fin, gzip.open(outfile, "wt") as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            
            seq = fin.readline()
            if not seq:
                sys.exit(f"Incomplete FASTA record: header without sequence in {infile}")
            
            seq = seq.rstrip("\n")
            base_header = header.rstrip("\n")

            # Find all TT positions
            positions = [i for i in range(len(seq) - 1) if seq[i : i + 2] == "TT"]
            seen_variants = set()

            # Generate one variant per TT position
            if positions:
                for idx, pos in enumerate(positions, start=1):
                    new_seq = seq[:pos] + "CC" + seq[pos + 2 :]
                    if new_seq in seen_variants:
                        continue
                    seen_variants.add(new_seq)
                    fout.write(f"{base_header} TT2CC_{idx}\n")
                    fout.write(new_seq + "\n")

            # Write the original sequence once at the end
            fout.write(base_header + "\n")
            fout.write(seq + "\n")


if __name__ == "__main__":
    main()
