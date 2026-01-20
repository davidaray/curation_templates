import argparse
import os
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description='Process genome files.')
    parser.add_argument('-g', '--good', required=True, help='Path to ${GENOME}_good.families.txt')
    parser.add_argument('-p', '--perfect', required=True, help='Path to ${GENOME}_perfect.families.txt')
    parser.add_argument('-c', '--cons', required=True, help='Path to ${GENOME}_confident_TE.cons.fa')
    parser.add_argument('-f', '--final', required=True, help='Path to ${GENOME}_file_final.0.1.txt')

    return parser.parse_args()

def load_file_to_dict(file_path):
    with open(file_path) as f:
        return {line.strip().split()[0]: line.strip().split() for line in f}

def main():
    args = parse_args()
    
    genome_name = os.path.basename(args.perfect).split('_')[0]
    
    good_families = load_file_to_dict(args.good)
    perfect_families = load_file_to_dict(args.perfect)
    file_final = load_file_to_dict(args.final)
    
    perfect_equivalencies = []
    good_equivalencies = []
    
    # Process perfect families
    with open(f"{genome_name}_perfect_equivalencies.txt", 'w') as pe_file:
        for key in perfect_families:
            if key in file_final:
                merge = f"{file_final[key][0]}#{file_final[key][1]}"
                perfect_equivalencies.append(f"{merge}\t{file_final[key][4]}")
                pe_file.write(f"{merge}\t{file_final[key][4]}\n")
    
    # Copy ${GENOME}_confident_TE.cons.fa to ${GENOME}_confident_TE.cons.trimmed.fa
    SeqIO.write(SeqIO.parse(args.cons, "fasta"), f"{genome_name}_confident_TE.cons.trimmed.fa", "fasta")
    
    trimmed_seqs = {record.id: record for record in SeqIO.parse(f"{genome_name}_confident_TE.cons.trimmed.fa", "fasta")}
    
    # Remove perfect matches from ${GENOME}_confident_TE.cons.trimmed.fa
    for key in perfect_families:
        if key in trimmed_seqs:
            del trimmed_seqs[key]
    
    SeqIO.write(trimmed_seqs.values(), f"{genome_name}_confident_TE.cons.trimmed.fa", "fasta")
    
    # Process good families and modify ${GENOME}_confident_TE.cons.trimmed.fa
    good_matches = []
    with open(f"{genome_name}_good_equivalencies.txt", 'w') as ge_file:
        for key in good_families:
            if key in file_final:
                merge = f"{file_final[key][0]}#{file_final[key][1]}"
                good_equivalencies.append(f"{merge}\t{file_final[key][4]}")
                ge_file.write(f"{merge}\t{file_final[key][4]}\n")
                if file_final[key][4] in trimmed_seqs:
                    record = trimmed_seqs.pop(file_final[key][4])
                    record.id = merge.replace('#', f"_{genome_name}#")
                    record.description = ''
                    good_matches.append(record)
                    print(f"Replaced header for {file_final[key][4]} with {record.id}")
    
    SeqIO.write(trimmed_seqs.values(), f"{genome_name}_confident_TE.cons.trimmed.fa", "fasta")
    SeqIO.write(good_matches, f"{genome_name}_good_matches.cons.fa", "fasta")
    
    print("Processing complete.")

if __name__ == "__main__":
    main()

