#!/usr/bin/env python3
"""
TE Header Modifier Script

This script modifies FASTA headers based on DeepTE classifications and a key file.
"""

import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Modify TE FASTA headers based on DeepTE calls and classification key'
    )
    parser.add_argument(
        '-if', '--input_fasta',
        required=True,
        help='Input FASTA file with TE consensus sequences'
    )
    parser.add_argument(
        '-mf', '--modified_fasta',
        required=True,
        help='Output FASTA file with modified headers'
    )
    parser.add_argument(
        '-uf', '--unmodified_fasta',
        required=True,
        help='Output FASTA file with unmodified headers'
    )
    parser.add_argument(
        '-af', '--all_fasta',
        required=True,
        help='Output FASTA file with all sequences (modified and unmodified)'
    )
    parser.add_argument(
        '-t', '--tsv',
        required=True,
        help='TSV file with original_header and deepte_call columns'
    )
    parser.add_argument(
        '-k', '--key',
        required=True,
        help='Key file mapping DeepTE calls to TE classifications'
    )
    
    return parser.parse_args()

def load_tsv(tsv_file):
    """Load TSV file with original headers and DeepTE calls"""
    deepte_calls = {}
    
    print(f"Loading DeepTE calls from: {tsv_file}")
    
    try:
        with open(tsv_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) >= 2:
                    original_header = parts[0]
                    deepte_call = parts[1]
                    deepte_calls[original_header] = deepte_call
        
        print(f"  Loaded {len(deepte_calls)} DeepTE calls")
        
    except FileNotFoundError:
        print(f"Error: TSV file '{tsv_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading TSV file: {e}")
        sys.exit(1)
    
    return deepte_calls

def load_key(key_file):
    """Load key file mapping DeepTE calls to TE classifications"""
    classification_key = {}
    
    print(f"Loading classification key from: {key_file}")
    
    try:
        with open(key_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t')
                if len(parts) >= 2:
                    deepte_class = parts[0]
                    te_class = parts[1]
                    classification_key[deepte_class] = te_class
        
        print(f"  Loaded {len(classification_key)} classification mappings")
        
    except FileNotFoundError:
        print(f"Error: Key file '{key_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading key file: {e}")
        sys.exit(1)
    
    return classification_key

def modify_header(original_header, deepte_call, classification_key):
    """
    Modify a header based on DeepTE call and classification key
    Returns (new_header, was_modified)
    """
    # Check if we should modify this header
    if deepte_call == "unknown" or deepte_call not in classification_key:
        return original_header, False
    
    # Get the new classification
    new_classification = classification_key[deepte_call]
    
    # Parse original header (format: ID#Classification)
    if '#' in original_header:
        seq_id, old_classification = original_header.split('#', 1)
        new_header = f"{seq_id}#{new_classification}"
        return new_header, True
    else:
        # If no # in header, just append the new classification
        new_header = f"{original_header}#{new_classification}"
        return new_header, True

def process_sequences(input_fasta, deepte_calls, classification_key):
    """
    Process all sequences and categorize them
    Returns: (modified_records, unmodified_records, all_records)
    """
    modified_records = []
    unmodified_records = []
    all_records = []
    
    print(f"\nProcessing sequences from: {input_fasta}")
    
    modified_count = 0
    unmodified_count = 0
    not_in_tsv = 0
    
    try:
        for record in SeqIO.parse(input_fasta, "fasta"):
            original_id = record.id
            original_description = record.description
            
            # Check if this sequence is in the DeepTE calls
            if original_id in deepte_calls:
                deepte_call = deepte_calls[original_id]
                
                # Try to modify the header
                new_id, was_modified = modify_header(original_id, deepte_call, classification_key)
                
                if was_modified:
                    # Create new record with modified header
                    new_record = SeqRecord(
                        record.seq,
                        id=new_id,
                        description=""
                    )
                    modified_records.append(new_record)
                    all_records.append(new_record)
                    modified_count += 1
                    print(f"  Modified: {original_id} -> {new_id}")
                else:
                    # Keep original record
                    unmodified_records.append(record)
                    all_records.append(record)
                    unmodified_count += 1
            else:
                # Not in TSV, keep as unmodified
                unmodified_records.append(record)
                all_records.append(record)
                not_in_tsv += 1
                unmodified_count += 1
    
    except FileNotFoundError:
        print(f"Error: Input FASTA file '{input_fasta}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing FASTA file: {e}")
        sys.exit(1)
    
    print(f"\nProcessing summary:")
    print(f"  Total sequences: {len(all_records)}")
    print(f"  Modified: {modified_count}")
    print(f"  Unmodified: {unmodified_count}")
    print(f"  Not in TSV: {not_in_tsv}")
    
    return modified_records, unmodified_records, all_records

def write_fasta(output_file, records):
    """Write sequences to FASTA file using BioPython"""
    try:
        count = SeqIO.write(records, output_file, "fasta")
        print(f"  Wrote {count} sequences to {output_file}")
    except Exception as e:
        print(f"Error writing FASTA file {output_file}: {e}")
        sys.exit(1)

def main():
    """Main function"""
    args = parse_arguments()
    
    # Load DeepTE calls and classification key
    deepte_calls = load_tsv(args.tsv)
    classification_key = load_key(args.key)
    
    # Process sequences
    modified_records, unmodified_records, all_records = process_sequences(
        args.input_fasta,
        deepte_calls,
        classification_key
    )
    
    # Write output files
    print("\nWriting output files:")
    write_fasta(args.modified_fasta, modified_records)
    write_fasta(args.unmodified_fasta, unmodified_records)
    write_fasta(args.all_fasta, all_records)
    
    print("\nProcessing complete!")

if __name__ == "__main__":
    main()
