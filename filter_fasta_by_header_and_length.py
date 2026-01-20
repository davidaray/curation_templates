import argparse
from Bio import SeqIO

def filter_fasta(input_file, query_string, max_length, removed_file, filtered_file):
    removed_sequences = []
    filtered_sequences = []

    # Process each sequence in the input FASTA file
    for record in SeqIO.parse(input_file, "fasta"):
        header = record.description
        sequence_length = len(record.seq)
        
        # Check if the sequence meets the criteria for removal
        if query_string in header and sequence_length > max_length:
            removed_sequences.append(record)
        else:
            filtered_sequences.append(record)

    # Save removed sequences
    with open(removed_file, "w") as rf:
        SeqIO.write(removed_sequences, rf, "fasta")

    # Save filtered sequences
    with open(filtered_file, "w") as ff:
        SeqIO.write(filtered_sequences, ff, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Filter sequences from a FASTA file.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file.")
    parser.add_argument("-q", "--query", required=True, help="String to search for in the header.")
    parser.add_argument("-l", "--length", type=int, required=True, help="Maximum length of sequences to retain.")
    parser.add_argument("-o", "--output", required=True, help="Base name for output files.")

    args = parser.parse_args()

    input_file = args.input
    query_string = args.query
    max_length = args.length
    base_output = args.output

    removed_file = f"{base_output}_removed.fasta"
    filtered_file = f"{base_output}_filtered.fasta"

    filter_fasta(input_file, query_string, max_length, removed_file, filtered_file)

    print(f"Filtered sequences saved to {filtered_file}")
    print(f"Removed sequences saved to {removed_file}")

if __name__ == "__main__":
    main()

