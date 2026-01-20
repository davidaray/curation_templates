import argparse
from Bio import SeqIO
from collections import defaultdict

def handle_duplicates(input_file, output_file):
    sequences = defaultdict(list)
    headers_seen = defaultdict(int)

    # Parse the input FASTA file
    for record in SeqIO.parse(input_file, "fasta"):
        # Normalize header by removing everything after the first occurrence of '#'
        header_base = record.id.split('#')[0]

        # Check if the sequence already exists in the sequences dictionary
        if str(record.seq) in sequences:
            sequences[str(record.seq)].append(record)
        else:
            # Handle identical headers with different sequences
            headers_seen[header_base] += 1
            if headers_seen[header_base] > 1:
                record.id = f"{header_base}-{headers_seen[header_base]}#{record.id.split('#', 1)[1]}"
                record.description = record.id

            sequences[str(record.seq)].append(record)

    # Write the non-duplicate sequences to the output file
    with open(output_file, "w") as output_handle:
        for seq_records in sequences.values():
            SeqIO.write(seq_records[0], output_handle, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Identify and remove duplicate sequences in a fasta file, renaming headers if needed.")
    parser.add_argument('-i', '--input_file', required=True, help="Path to the input fasta file.")
    parser.add_argument('-o', '--output_file', required=True, help="Path to the output fasta file without duplicates.")
    
    args = parser.parse_args()

    handle_duplicates(args.input_file, args.output_file)

if __name__ == "__main__":
    main()

