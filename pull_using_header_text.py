import argparse
from Bio import SeqIO
import fnmatch

def load_search_terms(search_file):
    """Loads search terms from a text file, one per line."""
    with open(search_file, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def filter_sequences(input_fasta, output_fasta, search_terms):
    """Filters sequences based on wildcard-supported search terms found anywhere in the header."""
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for record in SeqIO.parse(infile, 'fasta'):
            header = record.description
            # Check if any of the wildcard search terms match anywhere in the header
            if any(fnmatch.fnmatch(header, term) for term in search_terms):
                SeqIO.write(record, outfile, 'fasta')

def main():
    parser = argparse.ArgumentParser(description='Extract sequences with wildcard-supported text anywhere in the header.')
    parser.add_argument('-i', '--input', required=True, help='Input fasta file')
    parser.add_argument('-o', '--output', required=True, help='Output fasta file')
    parser.add_argument('-s', '--search_file', required=True, help='File containing wildcard search terms, one per line')

    args = parser.parse_args()

    # Load search terms from the provided file
    search_terms = load_search_terms(args.search_file)
    
    # Filter sequences based on the wildcard-supported search terms
    filter_sequences(args.input, args.output, search_terms)

if __name__ == "__main__":
    main()
