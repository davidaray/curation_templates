#!/usr/bin/env python3
"""
Modify transposable element consensus sequence headers based on classification results.

This script takes a FASTA file with TE consensus sequences and modifies headers
based on classification information provided in a TSV file and a key mapping file.
"""

import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Modify TE consensus sequence headers based on classification results'
    )
    parser.add_argument(
        '-if', '--input_fasta',
        required=True,
        help='Input FASTA file with TE consensus sequences'
    )
    parser.add_argument(
        '-mf', '--modified_fasta',
        required=True,
        help='Output FASTA file with modified headers only'
    )
    parser.add_argument(
        '-uf', '--unmodified_fasta',
        required=True,
        help='Output FASTA file with unmodified headers only'
    )
    parser.add_argument(
        '-af', '--all_fasta',
        required=True,
        help='Output FASTA file with all sequences (modified and unmodified)'
    )
    parser.add_argument(
        '-t', '--tsv',
        required=True,
        help='TSV file with three columns: known_header, TE_class, probability'
    )
    parser.add_argument(
        '-k', '--key',
        required=True,
        help='Key file mapping TE class to full classification (e.g., SINE -> SINE/SINE)'
    )
    
    return parser.parse_args()

def load_key_mapping(key_file):
    """
    Load the key file that maps TE class to full classification.
    
    Returns:
        dict: Mapping from TE class to full classification path
    """
    key_mapping = {}
    
    try:
        with open(key_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) >= 2:
                    te_class = parts[0].strip()
                    full_classification = parts[1].strip()
                    key_mapping[te_class] = full_classification
        
        print(f"Loaded {len(key_mapping)} TE class mappings from key file")
        return key_mapping
        
    except FileNotFoundError:
        print(f"Error: Key file '{key_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading key file: {e}")
        sys.exit(1)

def load_tsv_classifications(tsv_file):
    """
    Load the TSV file with classification results.
    
    Returns:
        dict: Mapping from known header to (TE_class, probability)
    """
    classifications = {}
    
    try:
        with open(tsv_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) >= 3:
                    known_header = parts[0].strip()
                    te_class = parts[1].strip()
                    probability = float(parts[2].strip())
                    classifications[known_header] = (te_class, probability)
        
        print(f"Loaded {len(classifications)} classifications from TSV file")
        return classifications
        
    except FileNotFoundError:
        print(f"Error: TSV file '{tsv_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading TSV file: {e}")
        sys.exit(1)

def modify_header(original_header, te_class, key_mapping):
    """
    Modify a header based on the TE class and key mapping.
    
    Args:
        original_header: Original header string (e.g., "amAme01.1.1#Unknown")
        te_class: TE class from classification (e.g., "SINE")
        key_mapping: Dictionary mapping TE class to full classification
    
    Returns:
        str: Modified header (e.g., "amAme01.1.1#SINE/SINE")
    """
    # Split header at '#' to get base name and classification
    if '#' in original_header:
        base_name, old_classification = original_header.split('#', 1)
        
        # Look up the new classification in the key mapping
        if te_class in key_mapping:
            new_classification = key_mapping[te_class]
            modified_header = f"{base_name}#{new_classification}"
            return modified_header
        else:
            print(f"Warning: TE class '{te_class}' not found in key mapping. Keeping original header.")
            return original_header
    else:
        print(f"Warning: Header '{original_header}' does not contain '#'. Keeping original.")
        return original_header

def process_sequences(input_fasta, classifications, key_mapping):
    """
    Process all sequences and determine which headers need modification.
    
    Returns:
        tuple: (modified_sequences, unmodified_sequences, all_sequences)
    """
    modified_sequences = []
    unmodified_sequences = []
    all_sequences = []
    
    modified_count = 0
    unmodified_count = 0
    
    try:
        for record in SeqIO.parse(input_fasta, "fasta"):
            original_header = record.description
            original_id = record.id
            
            # Check if this sequence needs modification
            if original_header in classifications:
                te_class, probability = classifications[original_header]
                
                # Modify the header
                new_header = modify_header(original_header, te_class, key_mapping)
                
                # Create new record with modified header
                new_record = SeqRecord(
                    record.seq,
                    id=new_header.split()[0],
                    description=new_header
                )
                
                modified_sequences.append(new_record)
                all_sequences.append(new_record)
                modified_count += 1
                
                print(f"Modified: {original_header} -> {new_header}")
                
            else:
                # Keep original header
                unmodified_sequences.append(record)
                all_sequences.append(record)
                unmodified_count += 1
        
        print(f"\nProcessing complete:")
        print(f"  Modified sequences: {modified_count}")
        print(f"  Unmodified sequences: {unmodified_count}")
        print(f"  Total sequences: {modified_count + unmodified_count}")
        
        return modified_sequences, unmodified_sequences, all_sequences
        
    except FileNotFoundError:
        print(f"Error: Input FASTA file '{input_fasta}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing sequences: {e}")
        sys.exit(1)

def write_fasta_file(sequences, output_file, description=""):
    """Write sequences to a FASTA file"""
    try:
        if sequences:
            SeqIO.write(sequences, output_file, "fasta")
            print(f"Wrote {len(sequences)} sequences to {output_file}")
        else:
            print(f"Warning: No sequences to write to {output_file}")
            # Create empty file
            open(output_file, 'w').close()
            
    except Exception as e:
        print(f"Error writing to {output_file}: {e}")
        sys.exit(1)

def main():
    """Main function"""
    args = parse_arguments()
    
    print("=" * 60)
    print("TE Header Modification Script")
    print("=" * 60)
    
    # Load key mapping
    print("\nStep 1: Loading key mapping...")
    key_mapping = load_key_mapping(args.key)
    
    # Load TSV classifications
    print("\nStep 2: Loading TSV classifications...")
    classifications = load_tsv_classifications(args.tsv)
    
    # Process sequences
    print("\nStep 3: Processing sequences...")
    modified_seqs, unmodified_seqs, all_seqs = process_sequences(
        args.input_fasta, 
        classifications, 
        key_mapping
    )
    
    # Write output files
    print("\nStep 4: Writing output files...")
    write_fasta_file(modified_seqs, args.modified_fasta, "modified sequences")
    write_fasta_file(unmodified_seqs, args.unmodified_fasta, "unmodified sequences")
    write_fasta_file(all_seqs, args.all_fasta, "all sequences")
    
    print("\n" + "=" * 60)
    print("Processing complete!")
    print("=" * 60)

if __name__ == "__main__":
    main()
