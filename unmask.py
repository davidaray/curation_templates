import argparse
import gzip
from Bio import SeqIO

def convert_to_uppercase(input_file, output_file):
    with gzip.open(input_file, 'rt') as handle_in:
        with gzip.open(output_file, 'wt') as handle_out:
            for record in SeqIO.parse(handle_in, 'fasta'):
                # Convert sequence to uppercase
                record.seq = record.seq.upper()
                # Write the updated sequence to the output file
                SeqIO.write(record, handle_out, 'fasta')

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Convert all nucleotide sequences in a .fa.gz file to uppercase.")
    parser.add_argument('-i', '--input', required=True, help="Input .fa.gz file containing nucleotide sequences")
    parser.add_argument('-o', '--output', required=True, help="Output .fa.gz file for the uppercase sequences")

    args = parser.parse_args()

    # Call the conversion function with provided arguments
    convert_to_uppercase(args.input, args.output)

if __name__ == "__main__":
    main()

