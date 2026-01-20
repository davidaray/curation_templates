import argparse
from Bio import SeqIO

def extract_sequences(input_fasta, terms, output_fasta):
    terms_set = set(terms)  # Convert list to set for faster lookups
    with open(output_fasta, 'w') as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            header_fields = record.description.split()  # Split header by whitespace
            if len(header_fields) > 1 and header_fields[1] in terms_set:
                SeqIO.write(record, out_f, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Extract sequences based on second field in FASTA headers.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-l", "--list", required=True, help="Comma-separated list of terms to match")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    args = parser.parse_args()

    terms = args.list.split(",")  # Convert list from string to Python list
    extract_sequences(args.input, terms, args.output)

if __name__ == "__main__":
    main()

