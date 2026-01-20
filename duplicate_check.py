import argparse
from Bio import SeqIO
from collections import defaultdict
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Process genome library files.')
    parser.add_argument('-i', '--input', required=True, help='Path to ${GENOME}_preliminary_library.fa')
    return parser.parse_args()

def find_duplicates_and_rename(input_file, genome_name):
    seq_dict = defaultdict(list)
    header_dict = defaultdict(list)
    duplicates = []
    renamed_records = []

    # Read sequences and group by sequence (check both orientations) and header
    for record in SeqIO.parse(input_file, "fasta"):
        seq = str(record.seq)
        rev_comp_seq = str(record.seq.reverse_complement())
        if seq in seq_dict or rev_comp_seq in seq_dict:
            duplicates.append(record)
        else:
            seq_dict[seq].append(record)
            seq_dict[rev_comp_seq].append(record)

        header_dict[record.id].append(record)

    # Write duplicates to file
    duplicates_file = f"{genome_name}_duplicates.fa"
    with open(duplicates_file, "w") as dup_handle:
        SeqIO.write(duplicates, dup_handle, "fasta")

    # Rename entries with duplicated headers
    for header, records in header_dict.items():
        if len(records) > 1:
            for i, record in enumerate(records, 1):
                new_header = f"{header.split('#')[0]}-{i}#{'#'.join(header.split('#')[1:])}"
                record.id = new_header
                record.description = ''
                renamed_records.append(record)
        else:
            renamed_records.append(records[0])

    # Write the final renamed library
    final_library_file = f"{genome_name}_final_library.fa"
    with open(final_library_file, "w") as final_handle:
        SeqIO.write(renamed_records, final_handle, "fasta")

    print(f"Processing complete. Duplicates written to {duplicates_file}. Final library written to {final_library_file}.")

def main():
    args = parse_args()
    input_file = args.input
    genome_name = os.path.basename(input_file).split('_')[0]
    find_duplicates_and_rename(input_file, genome_name)

if __name__ == "__main__":
    main()

"""
This Python script processes a genome library file in FASTA format, identifying duplicate sequences and renaming sequences with duplicated headers. Here's a detailed breakdown:

Key Functions
parse_args():

Uses argparse to parse the command-line argument -i or --input, which specifies the path to the input FASTA file.
Returns the parsed arguments.
find_duplicates_and_rename(input_file, genome_name):

Purpose: Processes the input FASTA file to identify and handle duplicate sequences and headers.
Steps:
Read sequences:
Reads sequences from the FASTA file.
Stores sequences in seq_dict (grouped by sequence and reverse complement).
Tracks headers in header_dict.
Identify duplicates:
Sequences that match (or reverse complement match) others in seq_dict are added to a duplicates list.
Rename entries with duplicated headers:
If a header occurs more than once in the file, appends a numeric suffix (-1, -2, etc.) to distinguish them.
Maintains a consistent naming convention for renamed headers.
Output results:
Writes duplicate sequences to a file named <genome_name>_duplicates.fa.
Writes the renamed library to a file named <genome_name>_final_library.fa.
Output Messages:
Prints confirmation messages showing the names of the output files.
main():

Parses the input arguments using parse_args().
Extracts the genome name (prefix of the input file name before the first _).
Calls find_duplicates_and_rename() with the input file and genome name.
Expected Behavior
Input:

A FASTA file containing sequences with potential duplicates or sequences with duplicate headers.
Example: GENOME_preliminary_library.fa
Output:

Duplicates file:
<genome_name>_duplicates.fa
Contains sequences that are identical (or reverse complement) to other sequences in the input file.
Final library file:
<genome_name>_final_library.fa
Contains all sequences with renamed headers for duplicated ones.
Printed message confirming the locations of the output files.
"""