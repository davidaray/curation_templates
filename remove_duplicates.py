import argparse
from Bio import SeqIO
from collections import defaultdict

def remove_duplicates(input_file, output_dir):
    # Dictionaries to store sequences and duplicates
    seq_dict = defaultdict(list)
    duplicates = []

    # Read the input FASTA file and populate the dictionary
    for record in SeqIO.parse(input_file, "fasta"):
        seq = str(record.seq)
        header = record.id

        if seq in seq_dict:
            # If the sequence is already in the dictionary, mark it as a duplicate
            duplicates.append(record)
        else:
            # Otherwise, store the record in the dictionary
            seq_dict[seq].append(record)

    # Write the unique sequences back to the input file
    with open(input_file, "w") as output_handle:
        for records in seq_dict.values():
            SeqIO.write(records, output_handle, "fasta")

    # Write the duplicates to the separate file
    output_id = input_file.split('/')[-1].split('_')[2]  # Extract ID from filename
    duplicates_file = f"{output_dir}/duplicated_sequences_{output_id}.pri.fa"
    with open(duplicates_file, "w") as duplicates_handle:
        SeqIO.write(duplicates, duplicates_handle, "fasta")

    print(f"Finished processing. Duplicates saved to {duplicates_file}.")

def main():
    parser = argparse.ArgumentParser(description="Remove exact duplicate sequences from a FASTA file.")
    parser.add_argument('-i', '--input', required=True, help="Path to the input FASTA file.")
    parser.add_argument('-o', '--output_dir', required=True, help="Directory to save the duplicate sequences file.")
    
    args = parser.parse_args()
    
    remove_duplicates(args.input, args.output_dir)

if __name__ == "__main__":
    main()

"""
This script removes exact duplicate sequences from a FASTA file and saves them to a separate file for review or further processing. Here's a breakdown of what it does:

Functionality Overview
Input FASTA File:

The script reads a FASTA file specified via the -i or --input argument.
It processes each sequence in the file, identifies exact duplicates (sequences with identical nucleotide content), and saves:
Unique sequences: back to the input file (overwrites it).
Duplicate sequences: to a separate file in a specified output directory.
Output Files:

Original Input File: Updated to retain only unique sequences.
Duplicates File: Contains all the duplicate sequences removed from the input file, named duplicated_sequences_<ID>.pri.fa, where <ID> is extracted from the input filename.
Detailed Breakdown
remove_duplicates(input_file, output_dir):

Step 1: Parse Input File:

Reads all sequences in the input FASTA file using Biopython's SeqIO.
Stores unique sequences in a dictionary seq_dict (keyed by sequence).
Tracks duplicate sequences in a list duplicates.
Step 2: Identify Duplicates:

If a sequence is already in seq_dict, it is added to the duplicates list.
If a sequence is not in seq_dict, it is added as a unique sequence.
Step 3: Write Outputs:

Writes the unique sequences back to the input file, overwriting it.
Writes the duplicate sequences to a new file in the specified output directory.
Filename Parsing:

Extracts an identifier (<ID>) from the input filename by splitting it on underscores (_) and taking the third segment (input_file.split('/')[-1].split('_')[2]).
Example: From path/to/file_GENOME_123_library.fa, it would extract 123.
Output Files:

Unique Sequences: Overwrites the input FASTA file.
Duplicates: Saved as duplicated_sequences_<ID>.pri.fa in the specified output directory.
main():

Uses argparse to parse command-line arguments:
-i / --input: Path to the input FASTA file.
-o / --output_dir: Directory where the duplicate sequences file will be saved.
Calls remove_duplicates() with the parsed arguments.
"""