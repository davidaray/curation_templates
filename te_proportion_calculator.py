import gzip
import argparse
import os

def get_genome_size(summary_file):
    with gzip.open(summary_file, 'rt') as f:
        for line in f:
            if "Total Length:" in line:
                parts = line.split()
                for part in parts:
                    if part.isdigit():
                        return int(part)
    return 0

def parse_repeatmasker(out_file, threshold):
    te_bp = {"LINE": 0, "SINE": 0, "LTR": 0, "DNA": 0, "RC": 0}
    with gzip.open(out_file, 'rt') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue  # Skip headers and empty lines
            fields = line.split()
            if len(fields) < 15:
                continue  # Ensure correct format
            try:
                divergence = float(fields[1])
                te_class = fields[10].split('/')[0]  # Extract TE class
                te_length = int(fields[6]) - int(fields[5])  # Calculate bp occupied
                if divergence < threshold and te_class in te_bp:
                    te_bp[te_class] += te_length  # Add bp occupied to TE class
            except ValueError:
                continue  # Skip lines with unexpected formats
    return te_bp

def process_files(directory, mapping_file, threshold, output_file):
    basename = os.path.splitext(output_file)[0]
    output_file_extended = basename + "_extended.tsv" 
    mapping = {}
    with open(mapping_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            fields = line.strip().split('\t')
            mapping[fields[2]] = fields[1]  # Species_ID to Binomial_Species_Name
    
    results = []
    extended_results = []
    for species_id in mapping:
        out_file = os.path.join(directory, f"{species_id}.fa.out.gz")
        summary_file = os.path.join(directory, f"{species_id}.summary.gz")
        print('Processing: ', species_id)
        
        if not os.path.exists(out_file):
            print(out_file, "does not exist")
            continue
        if not os.path.exists(summary_file):
            print(summary_file, "does not exist")
            continue
        
        genome_size = get_genome_size(summary_file)
        if genome_size == 0:
            print("Could not determine genome size for", species_id)
            continue
        
        te_totals = parse_repeatmasker(out_file, threshold)
        for te_class, total_bp in te_totals.items():
            proportion = total_bp / genome_size if genome_size > 0 else 0
            results.append([mapping[species_id], te_class, proportion])
            extended_results.append([mapping[species_id], te_class, proportion, total_bp, genome_size])
        
    with open(output_file, 'w') as f:
        f.write("SPECIES\tTE\tPROP\n")
        for result in results:
            f.write("\t".join(map(str, result)) + "\n")
            
    with open(output_file_extended, 'w') as f:
        f.write("SPECIES\tTE\tPROP\tTOTALBP\tGENOMESIZE\n")
        for extended_result in extended_results:
            f.write("\t".join(map(str, extended_result)) + "\n")
            
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", required=True, help="Directory containing input files")
    parser.add_argument("-m", "--mapping_file", required=True, help="Mapping file")
    parser.add_argument("-t", "--threshold", type=float, required=True, help="Maximum divergence threshold")
    parser.add_argument("-o", "--output_file", required=True, help="Output TSV file")
    args = parser.parse_args()
    
    process_files(args.directory, args.mapping_file, args.threshold, args.output_file)

if __name__ == "__main__":
    main()
    
