import argparse
import csv

def parse_fasta(fasta_file):
    """Parse the headers and sequences of a FASTA file into a dictionary."""
    headers = {}
    sequences = {}
    with open(fasta_file, 'r') as f:
        header = None
        sequence_lines = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    headers[header] = header_info
                    sequences[header] = ''.join(sequence_lines)
                header_full = line[1:]  # Remove '>'
                header, header_info = header_full.split('#', 1)
                sequence_lines = []
            else:
                sequence_lines.append(line)
        # Add last sequence
        if header:
            headers[header] = header_info
            sequences[header] = ''.join(sequence_lines)
    return headers, sequences

def compare_headers(deepte_headers, teclass2_headers, deepte_sequences, mismatch_file, match_file, output_fasta):
    """Compare headers between DeepTE and TEclass2, and write mismatches, matches, and sequences to files."""
    
    with open(mismatch_file, 'w', newline='') as mismatch_csv, open(match_file, 'w', newline='') as match_csv, open(output_fasta, 'w') as fasta_out:
        fieldnames = ['header id', 'DeepTE header', 'TEclass2 header']
        
        # Create CSV writers for mismatches and matches
        mismatch_writer = csv.DictWriter(mismatch_csv, fieldnames=fieldnames, delimiter='\t')
        match_writer = csv.DictWriter(match_csv, fieldnames=fieldnames, delimiter='\t')
        
        # Write headers to both output files
        mismatch_writer.writeheader()
        match_writer.writeheader()
        
        for header_id in deepte_headers:
            deepte_header = deepte_headers[header_id]
            teclass2_header = teclass2_headers.get(header_id)
            
            if teclass2_header:
                # Check for matches and mismatches
                if deepte_header == teclass2_header:
                    # Write match to match file
                    match_writer.writerow({
                        'header id': header_id,
                        'DeepTE header': deepte_header,
                        'TEclass2 header': teclass2_header
                    })
                    # Copy sequence to output FASTA with no changes
                    fasta_out.write(f'>{header_id}#{deepte_header}\n{deepte_sequences[header_id]}\n')
                else:
                    # Write mismatch to mismatch file
                    mismatch_writer.writerow({
                        'header id': header_id,
                        'DeepTE header': deepte_header,
                        'TEclass2 header': teclass2_header
                    })
                    # Modify header in output FASTA
                    fasta_out.write(f'>{header_id}#Unknown/Unknown\n{deepte_sequences[header_id]}\n')

def main():
    # Set up argparse to handle command line arguments
    parser = argparse.ArgumentParser(description="Compare headers between DeepTE and TEclass2 FASTA files and output sequences.")
    parser.add_argument('-d', '--deepte', required=True, help="Input DeepTE FASTA file")
    parser.add_argument('-t', '--teclass2', required=True, help="Input TEclass2 FASTA file")
    parser.add_argument('-m', '--mismatch_output', default='mismatches.tsv', help="Output file for mismatches (default: mismatches.tsv)")
    parser.add_argument('-p', '--match_output', default='matches.tsv', help="Output file for matches (default: matches.tsv)")
    parser.add_argument('-of', '--output_fasta', required=True, help="Output FASTA file for sequences")
    
    args = parser.parse_args()
    
    # Parse both FASTA files to get headers and sequences
    deepte_headers, deepte_sequences = parse_fasta(args.deepte)
    teclass2_headers, _ = parse_fasta(args.teclass2)
    
    # Compare headers and write both matches, mismatches, and sequences to output FASTA
    compare_headers(deepte_headers, teclass2_headers, deepte_sequences, args.mismatch_output, args.match_output, args.output_fasta)

if __name__ == '__main__':
    main()

