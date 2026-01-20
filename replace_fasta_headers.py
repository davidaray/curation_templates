#!/usr/bin/env python

import argparse
from Bio import SeqIO

def replace_headers(fasta_file, output_file, new_header):
    # Open the input FASTA file for reading
    with open(fasta_file, "r") as input_handle:
        # Open the output file for writing
        with open(output_file, "w") as output_handle:
            # Iterate over each sequence record in the input FASTA file
            for i, record in enumerate(SeqIO.parse(input_handle, "fasta"), 1):
                # Replace the header with the new header string and append a unique number
                record.id = f"{new_header}"
                record.description = ""  # Optional: Clear the description part
                # Write the modified sequence to the output file
                SeqIO.write(record, output_handle, "fasta")

def main():
    # Create an argument parser
    parser = argparse.ArgumentParser(description="Replace FASTA headers with a given string")

    # Define arguments
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file with modified headers")
    parser.add_argument("-n", "--new_header", required=True, help="New header string to replace existing headers")

    # Parse arguments
    args = parser.parse_args()

    # Call the function to replace headers
    replace_headers(args.input, args.output, args.new_header)

if __name__ == "__main__":
    main()

