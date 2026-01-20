#!/usr/bin/env python

from Bio.SeqUtils import six_frame_translations
from Bio.Data.CodonTable import TranslationError
from Bio import SeqIO
import argparse
import os

def get_sequence_from_file(file_path):
    """Reads sequences from the given file and returns a list of SeqRecord objects."""
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(record)
    return sequences

def process_sequences(sequences, output_dir, output_file):
    """Processes sequences, skips ones with TranslationError, and saves valid sequences."""
    valid_sequences = []
    
    for record in sequences:
        try:
            # Attempt translation of the sequence to catch translation errors
            _ = six_frame_translations(record.seq)
            valid_sequences.append(record)  # Save if no error occurred
        except TranslationError as e:
            print(f"Error in sequence {record.id}: {e}")
            # Log the error to a file
            with open(os.path.join(output_dir, "translation_errors.log"), "a") as log_file:
                log_file.write(f"Error in sequence {record.id}: {e}\n")
    
    # Save valid sequences to output file
    with open(output_file, "w") as output_handle:
        SeqIO.write(valid_sequences, output_handle, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Clean and process sequences, removing those with translation errors.")
    parser.add_argument("-i", dest="input_file", required=True, help="Input FASTA file containing sequences.")
    parser.add_argument("-o", dest="output_file", required=True, help="Output FASTA file to save valid sequences.")
    parser.add_argument("-d", dest="output_dir", default="./", help="Directory to store logs. Default: ./")
    args = parser.parse_args()

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # Read sequences from the input file
    sequences = get_sequence_from_file(args.input_file)

    # Process sequences and save valid ones
    process_sequences(sequences, args.output_dir, args.output_file)

if __name__ == "__main__":
    main()
