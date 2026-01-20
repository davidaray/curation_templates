import argparse
from Bio import SeqIO

def extract_sequences(input_fasta, query_id, output_fasta):
    """
    Extracts sequences from the input FASTA file whose ID exactly matches the query ID
    and writes them to an output FASTA file.

    Parameters:
    - input_fasta: Path to the input FASTA file
    - query_id: Exact ID to search for
    - output_fasta: Path to the output FASTA file
    """
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        count = 0
        for record in SeqIO.parse(infile, 'fasta'):
            if record.id == query_id:  # exact match
                SeqIO.write(record, outfile, 'fasta')
                count += 1
        print(f"Extracted {count} sequences with ID '{query_id}' into {output_fasta}.")

def main():
    parser = argparse.ArgumentParser(
        description="Extract sequences from a FASTA file based on an exact ID match."
    )
    parser.add_argument('-i', '--input', required=True, help="Path to the input FASTA file.")
    parser.add_argument('-q', '--query', required=True, help="Exact sequence ID to search for.")
    parser.add_argument('-o', '--output', required=True, help="Path to the output FASTA file.")
    
    args = parser.parse_args()
    
    extract_sequences(args.input, args.query, args.output)

if __name__ == '__main__':
    main()
