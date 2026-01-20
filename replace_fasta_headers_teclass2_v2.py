#!/usr/bin/env python3
import argparse
import pandas as pd
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Replace FASTA headers using mapping from table")
    parser.add_argument("-i", "--input_table", required=True, help="Input TSV file with Name and SuperClass")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file with modified headers")
    args = parser.parse_args()

    # Load table
    mapping_df = pd.read_csv(args.input_table, sep="\t")
    name_to_class = dict(zip(mapping_df["Name"], mapping_df["SuperClass"]))

    # Process fasta
    records = []
    for record in SeqIO.parse(args.fasta, "fasta"):
        base_id = record.id.split("#")[0]
        if base_id in name_to_class:
            record.id = f"{base_id}#{name_to_class[base_id]}"
            record.description = ""  # clear to avoid duplication in FASTA
        records.append(record)

    # Write output
    SeqIO.write(records, args.output, "fasta")

if __name__ == "__main__":
    main()

