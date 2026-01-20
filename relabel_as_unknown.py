#!/usr/bin/env python3
"""
FASTA Unknown Labeler Script

This script replaces everything after '#' in FASTA headers with 'Unknown/Unknown'.
"""

import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Replace classification after # with Unknown/Unknown in FASTA headers'
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input FASTA file'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output FASTA file with modified headers'
    )
    
    return parser.parse_args()

def relabel_headers(input_fasta, output_fasta):
    """
    Process FASTA file and replace everything after # with Unknown/Unknown
    """
    print(f"Processing FASTA file: {input_fasta}")
    
    modified_records = []
    modified_count = 0
    no_hash_count = 0
    total_count = 0
    
    try:
        for record in SeqIO.parse(input_fasta, "fasta"):
            total_count += 1
            
            # Check if header contains '#'
            if '#' in record.id:
                # Split at '#' and replace everything after with Unknown/Unknown
                seq_id = record.id.split('#')[0]
                new_id = f"{seq_id}#Unknown/Unknown"
                
                # Create new record with modified header
                new_record = SeqRecord(
                    record.seq,
                    id=new_id,
                    description=""
                )
                
                modified_records.append(new_record)
                modified_count += 1
                
                if modified_count <= 5:  # Show first 5 examples
                    print(f"  {record.id} -> {new_id}")
            else:
                # No '#' in header, keep as is
                modified_records.append(record)
                no_hash_count += 1
        
        # Write output file
        print(f"\nWriting output to: {output_fasta}")
        SeqIO.write(modified_records, output_fasta, "fasta")
        
        # Print summary
        print(f"\nProcessing complete!")
        print(f"  Total sequences: {total_count}")
        print(f"  Modified (had #): {modified_count}")
        print(f"  Unchanged (no #): {no_hash_count}")
        
    except FileNotFoundError:
        print(f"Error: Input file '{input_fasta}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)

def main():
    """Main function"""
    args = parse_arguments()
    relabel_headers(args.input, args.output)

if __name__ == "__main__":
    main()
