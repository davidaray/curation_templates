import argparse
from Bio import SeqIO

def filter_sequences(input_fasta, output_fasta, search_term):
    """Filters sequences based on the specified text after the '#' in the header."""
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for record in SeqIO.parse(infile, 'fasta'):
            header = record.description
            # Check if '#' exists in the header and search term appears after it
            if '#' in header:
                header_after_hash = header.split('#', 1)[1]
                if search_term in header_after_hash:
                    SeqIO.write(record, outfile, 'fasta')

def main():
    parser = argparse.ArgumentParser(description='Extract sequences with a specified string after # in the header.')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output FASTA file')
    parser.add_argument('-s', '--search', required=True, help='Text string to search for after # in the header')

    args = parser.parse_args()

    # Filter sequences based on the search term
    filter_sequences(args.input, args.output, args.search)

if __name__ == "__main__":
    main()
