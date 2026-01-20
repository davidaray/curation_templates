import os
import argparse
import gzip
from pathlib import Path

def parse_faidx(file_path, num_scaffolds):
    scaffolds = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            scaffold_name = parts[0]
            scaffold_size = int(parts[1])
            scaffolds.append((scaffold_name, scaffold_size))
    # Sort scaffolds by size in descending order and select the top N scaffolds
    top_scaffolds = sorted(scaffolds, key=lambda x: x[1], reverse=True)[:num_scaffolds]
    return [scaffold[0] for scaffold in top_scaffolds]

def extract_scaffold_to_fa(input_file, scaffold_name, output_file):
    with gzip.open(input_file, 'rt') as infile, open(output_file, 'w') as outfile:
        capture = False
        for line in infile:
            if line.startswith(">"):
                if capture:
                    break
                capture = scaffold_name in line
            if capture:
                outfile.write(line)

def main(args):
    species_pairs = [line.strip() for line in open(args.list).readlines()]

    # Create the output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    
    for pair in species_pairs:
        hap1, hap2, num_scaffolds = pair.split()
        num_scaffolds = int(num_scaffolds)
        
        # Locate .fa.gz and faidx files for each haplotype
        hap1_fa = Path(args.directory) / f"{hap1}.fa.gz"
        hap2_fa = Path(args.directory) / f"{hap2}.fa.gz"
        hap1_faidx = Path(args.directory) / f"{hap1}.faidx"
        hap2_faidx = Path(args.directory) / f"{hap2}.faidx"
        
        # Parse faidx files to get the top N scaffolds by size
        top_scaffolds_hap1 = parse_faidx(hap1_faidx, num_scaffolds)
        top_scaffolds_hap2 = parse_faidx(hap2_faidx, num_scaffolds)
        
        # Extract corresponding scaffolds and create output .fa files
        for i, scaffold in enumerate(top_scaffolds_hap1):
            hap1_output = Path(args.output) / f"{hap1}.{scaffold}.fa"
            hap2_output = Path(args.output) / f"{hap2}.{top_scaffolds_hap2[i]}.fa"
            
            extract_scaffold_to_fa(hap1_fa, scaffold, hap1_output)
            extract_scaffold_to_fa(hap2_fa, top_scaffolds_hap2[i], hap2_output)
            print(f"Created {hap1_output} and {hap2_output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract top N largest scaffolds for paired haplotypes.")
    parser.add_argument("-l", "--list", required=True, help="File with list of paired species names and scaffold count.")
    parser.add_argument("-d", "--directory", required=True, help="Directory containing .fa.gz and .faidx files.")
    parser.add_argument("-o", "--output", required=True, help="Output directory for new .fa files.")
    
    args = parser.parse_args()
    main(args)
