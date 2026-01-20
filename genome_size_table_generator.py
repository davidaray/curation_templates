#!/usr/bin/env python3
"""
Genome Size Table Generator

This script reads a mapping file and genome assembly directory to create
a genome size table with species information.
"""

import argparse
import pandas as pd
import gzip
import os
import sys

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Generate genome size table from mapping file and genome assemblies'
    )
    parser.add_argument(
        '-m', '--mapping',
        required=True,
        help='Path to mapping TSV file with columns: Taxonomic_Family, Binomial_Species_Name, Species_ID, Assembly_ID, Old_ID, TOGA_bed'
    )
    parser.add_argument(
        '-d', '--directory',
        required=True,
        help='Directory containing genome assemblies named <Species_ID>.fa.gz'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output file name for genome size table'
    )
    
    return parser.parse_args()

def read_mapping_file(mapping_file):
    """Read the mapping TSV file and return as DataFrame"""
    try:
        print(f"Reading mapping file: {mapping_file}")
        df = pd.read_csv(mapping_file, sep='\t')
        
        # Check for required columns
        required_cols = ['Species_ID', 'Binomial_Species_Name']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            print(f"Error: Missing required columns in mapping file: {missing_cols}")
            sys.exit(1)
        
        print(f"Loaded {len(df)} entries from mapping file")
        return df
        
    except FileNotFoundError:
        print(f"Error: Mapping file '{mapping_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading mapping file: {e}")
        sys.exit(1)

def calculate_genome_size(fasta_file):
    """Calculate total genome size from a gzipped FASTA file"""
    total_size = 0
    
    try:
        with gzip.open(fasta_file, 'rt') as f:
            for line in f:
                line = line.strip()
                # Skip header lines (starting with '>')
                if not line.startswith('>'):
                    total_size += len(line)
        
        return total_size
        
    except FileNotFoundError:
        print(f"  Warning: File not found: {fasta_file}")
        return None
    except gzip.BadGzipFile:
        print(f"  Warning: Not a valid gzip file: {fasta_file}")
        return None
    except Exception as e:
        print(f"  Warning: Error reading {fasta_file}: {e}")
        return None

def generate_genome_size_table(mapping_df, genome_dir, output_file):
    """Generate genome size table from mapping and genome files"""
    
    results = []
    mutation_rate = "2.2E-09"  # Fixed mutation rate
    
    print(f"\nProcessing genome assemblies from: {genome_dir}")
    
    for idx, row in mapping_df.iterrows():
        species_id = row['Species_ID']
        binomial_name = row['Binomial_Species_Name']
        
        # Construct genome file path
        genome_file = os.path.join(genome_dir, f"{species_id}.fa.gz")
        
        print(f"Processing {species_id}...", end=' ')
        
        # Calculate genome size
        genome_size = calculate_genome_size(genome_file)
        
        if genome_size is not None:
            print(f"Genome size: {genome_size:,} bp")
            results.append({
                'Species_ID': species_id,
                'Genome_Size_bp': genome_size,
                'Mutation_Rate': mutation_rate,
                'Binomial_Species_Name': binomial_name
            })
        else:
            print(f"Skipped (file not found or error)")
    
    if not results:
        print("\nError: No genome sizes were successfully calculated.")
        sys.exit(1)
    
    # Create DataFrame and save
    results_df = pd.DataFrame(results)
    
    # Save without header and without index, tab-separated
    results_df.to_csv(output_file, sep='\t', index=False, header=False)
    
    print(f"\nGenome size table saved to: {output_file}")
    print(f"Successfully processed {len(results)} genomes")
    
    # Print summary statistics
    print("\nSummary:")
    print(f"  Minimum genome size: {results_df['Genome_Size_bp'].min():,} bp")
    print(f"  Maximum genome size: {results_df['Genome_Size_bp'].max():,} bp")
    print(f"  Average genome size: {results_df['Genome_Size_bp'].mean():,.0f} bp")

def main():
    """Main function"""
    args = parse_arguments()
    
    # Check if directory exists
    if not os.path.isdir(args.directory):
        print(f"Error: Directory '{args.directory}' does not exist.")
        sys.exit(1)
    
    # Read mapping file
    mapping_df = read_mapping_file(args.mapping)
    
    # Generate genome size table
    generate_genome_size_table(mapping_df, args.directory, args.output)
    
    print("\nDone!")

if __name__ == "__main__":
    main()
