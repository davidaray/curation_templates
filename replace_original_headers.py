import argparse
from Bio import SeqIO

def parse_arguments():
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Alter fasta headers using another fasta file.")
    parser.add_argument("-f1", "--first_fasta", required=True, help="First input fasta file with existing headers")
    parser.add_argument("-f2", "--second_fasta", required=True, help="Second input fasta file with altered headers")
    parser.add_argument("-o", "--output_fasta", help="Optional output file name")
    return parser.parse_args()

def extract_text1(header):
    """Extracts the text1 portion from a header of the format 'text1#text2/text3'."""
    return header.split('#')[0]

def load_altered_headers(second_fasta):
    """Loads altered headers from the second input fasta."""
    altered_headers = {}
    for record in SeqIO.parse(second_fasta, "fasta"):
        text1 = extract_text1(record.id)  # Extract 'text1' from the second file's header
        altered_headers[text1] = record.id  # Map text1 to the altered header
    return altered_headers

def alter_fasta_headers(first_fasta, altered_headers, output_fasta):
    """Alters headers in the first fasta file using the altered headers from the second file."""
    with open(output_fasta, "w") as out_f:
        for record in SeqIO.parse(first_fasta, "fasta"):
            text1 = extract_text1(record.id)
            # If text1 is in the altered headers, replace the header
            if text1 in altered_headers:
                record.id = altered_headers[text1]
                record.description = ""  # Clear the description field
            SeqIO.write(record, out_f, "fasta")

def main():
    args = parse_arguments()

    # Load altered headers from the second fasta file
    altered_headers = load_altered_headers(args.second_fasta)

    # Define output file name if not provided
    if not args.output_fasta:
        output_fasta = args.first_fasta.rsplit('.fa', 1)[0] + "_deep_te.fa"
    else:
        output_fasta = args.output_fasta

    # Alter headers in the first fasta file using the loaded altered headers
    alter_fasta_headers(args.first_fasta, altered_headers, output_fasta)

if __name__ == "__main__":
    main()

