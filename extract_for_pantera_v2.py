import os
import argparse
import gzip
import subprocess
import logging
from pathlib import Path

def setup_logging():
    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.INFO,
        handlers=[logging.StreamHandler()]
    )

def parse_faidx(file_path, min_size):
    scaffolds = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            scaffold_name = parts[0]
            scaffold_size = int(parts[1])
            if scaffold_size >= min_size:
                scaffolds.append((scaffold_name, scaffold_size))
                logging.info(f"Accepted scaffold '{scaffold_name}' with size {scaffold_size} bp (min size: {min_size} bp).")
            else:
                logging.info(f"Rejected scaffold '{scaffold_name}' with size {scaffold_size} bp (below min size threshold).")
    return sorted(scaffolds, key=lambda x: x[1], reverse=True)

def align_scaffolds(scaffold1_path, scaffold2_path, output_path, threads):
    command = f"minimap2 -x asm5 -t {threads} {scaffold1_path} {scaffold2_path} > {output_path}"
    logging.info(f"Running alignment: {command}")
    subprocess.run(command, shell=True, check=True)
    logging.info(f"Alignment complete. Output saved to {output_path}.")

def calculate_similarity(paf_file):
    total_matches, total_length = 0, 0
    with open(paf_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            match_len = int(parts[9])
            total_len = int(parts[10])
            total_matches += match_len
            total_length += total_len
    similarity = (total_matches / total_length) if total_length > 0 else 0
    logging.info(f"Calculated similarity: {similarity * 100:.2f}% from {paf_file}")
    return similarity

def extract_scaffold(input_file, scaffold_name, output_file):
    with gzip.open(input_file, 'rt') as infile, open(output_file, 'w') as outfile:
        capture = False
        for line in infile:
            if line.startswith(">"):
                if capture:
                    break
                capture = scaffold_name in line
            if capture:
                outfile.write(line)
    logging.info(f"Extracted scaffold '{scaffold_name}' from {input_file} to {output_file}")

def main(args):
    setup_logging()
    logging.info("Starting script to identify and extract paired autosomes from phased haplotypes.")

    species_pairs = [line.strip().split() for line in open(args.list).readlines()]
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    results = []

    for hap1, hap2 in species_pairs:
        logging.info(f"Processing species pair: {hap1} and {hap2}")

        hap1_fa = Path(args.directory) / f"{hap1}.fa.gz"
        hap2_fa = Path(args.directory) / f"{hap2}.fa.gz"
        hap1_faidx = Path(args.directory) / f"{hap1}.faidx"
        hap2_faidx = Path(args.directory) / f"{hap2}.faidx"

        hap1_scaffolds = parse_faidx(hap1_faidx, args.minsize)
        hap2_scaffolds = parse_faidx(hap2_faidx, args.minsize)

        for scaffold1, size1 in hap1_scaffolds:
            for scaffold2, size2 in hap2_scaffolds:
                size_diff = abs(size1 - size2) / max(size1, size2)
                logging.info(f"Comparing {scaffold1} (size: {size1}) and {scaffold2} (size: {size2}) - size differential: {size_diff:.4f}")
                
                if size_diff < args.sizediff:
                    logging.info(f"Size differential {size_diff:.4f} is below threshold {args.sizediff}. Proceeding with alignment.")
                    
                    scaffold1_fa = output_dir / f"{hap1}.{scaffold1}.fa"
                    scaffold2_fa = output_dir / f"{hap2}.{scaffold2}.fa"
                    extract_scaffold(hap1_fa, scaffold1, scaffold1_fa)
                    extract_scaffold(hap2_fa, scaffold2, scaffold2_fa)

                    alignment_output = output_dir / f"{hap1}_{scaffold1}_vs_{hap2}_{scaffold2}.paf"
                    align_scaffolds(scaffold1_fa, scaffold2_fa, alignment_output, args.threads)

                    similarity = calculate_similarity(alignment_output)
                    if similarity >= args.similarity:
                        logging.info(f"Similarity {similarity * 100:.2f}% meets threshold {args.similarity * 100:.2f}%.")
                        
                        paired_output = output_dir / f"{hap1}_{scaffold1}_{hap2}_{scaffold2}_paired.fa"
                        with open(paired_output, 'w') as outfile:
                            for file_path in [scaffold1_fa, scaffold2_fa]:
                                with open(file_path, 'r') as f:
                                    outfile.write(f.read())
                        logging.info(f"Paired autosomes saved to {paired_output}")

                        results.append((f"{hap1}_{scaffold1}_vs_{hap2}_{scaffold2}", size_diff, similarity * 100))
                    
                    scaffold1_fa.unlink(missing_ok=True)
                    scaffold2_fa.unlink(missing_ok=True)
                    alignment_output.unlink(missing_ok=True)
                else:
                    logging.info(f"Size differential {size_diff:.4f} exceeds threshold {args.sizediff}. Skipping alignment.")

    with open(output_dir / "homologous_scaffolds_summary.txt", 'w') as f:
        f.write("Filename\tSize Differential\tSimilarity (%)\n")
        for filename, size_diff, similarity in results:
            f.write(f"{filename}\t{size_diff:.4f}\t{similarity:.2f}\n")
    logging.info("Script completed. Summary of homologous scaffolds saved.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identify and extract paired autosomes from phased haplotypes.")
    parser.add_argument("-l", "--list", required=True, help="File with list of paired species names.")
    parser.add_argument("-d", "--directory", required=True, help="Directory containing .fa.gz and .faidx files.")
    parser.add_argument("-o", "--output", required=True, help="Output directory for paired autosome files.")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads for minimap2 alignment (default: 1).")
    parser.add_argument("--minsize", type=int, default=1000000, help="Minimum scaffold size to consider (default: 1,000,000 bp).")
    parser.add_argument("-s", "--similarity", type=float, default=0.9, help="Minimum similarity threshold for paired autosomes (default: 0.9).")
    parser.add_argument("--sizediff", type=float, default=0.05, help="Maximum size difference ratio for paired autosomes (default: 0.05).")

    args = parser.parse_args()
    main(args)
