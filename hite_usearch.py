import argparse
import os
import pandas as pd
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description='Process genome files.')
    parser.add_argument('-u', '--usearch', required=True, help='Path to ${GENOME}_usearch_80_hits.tsv')
    parser.add_argument('-c', '--cons', required=True, help='Path to ${GENOME}_confident_TE.cons.fa')
    parser.add_argument('-m', '--mammals', required=True, help='Path to mammals.plus.covid_bats2.04072022.fa')
    return parser.parse_args()

def process_usearch_file(usearch_file):
    df = pd.read_csv(usearch_file, sep='\t', header=None)
    df.columns = ['A', 'B', 'C', 'D', 'E']
    df['F'] = df['D'] / df['E']
    df['G'] = df.apply(lambda row: 'same' if row['C'] > 95 and row['D'] > 80 and 0.8 <= row['F'] <= 1.2 
                       else ('family' if 80 <= row['C'] <= 95 and row['D'] > 80 and 0.8 <= row['F'] <= 1.2 else ''), axis=1)
    df.to_csv(usearch_file, sep='\t', index=False, header=False)
    return df

def copy_and_trim_fasta(cons_file, trimmed_fasta):
    SeqIO.write(SeqIO.parse(cons_file, "fasta"), trimmed_fasta, "fasta")

def process_fasta(df, cons_file, same_file, family_file, genome_name):
    trimmed_seqs = {record.id: record for record in SeqIO.parse(cons_file, "fasta")}
    
    same_seqs = []
    family_seqs = []
    
    for _, row in df.iterrows():
        if row['G'] == 'same':
            if row['A'] in trimmed_seqs:
                same_seqs.append(trimmed_seqs.pop(row['A']))
        elif row['G'] == 'family':
            if row['A'] in trimmed_seqs:
                record = trimmed_seqs.pop(row['A'])
                new_header = row['B'].replace('#', f"_{genome_name}#")
                record.id = new_header
                record.description = ''
                family_seqs.append(record)
                print(f"Replaced header for {row['A']} with {record.id}")
    
    SeqIO.write(same_seqs, same_file, "fasta")
    SeqIO.write(family_seqs, family_file, "fasta")
    SeqIO.write(trimmed_seqs.values(), cons_file, "fasta")

def main():
    args = parse_args()
    
    genome_name = os.path.basename(args.usearch).split('_')[0]
    
    df = process_usearch_file(args.usearch)
    
    trimmed_fasta = f"{genome_name}_confident_TE.cons.trimmed.fa"
    copy_and_trim_fasta(args.cons, trimmed_fasta)
    
    same_file = f"{genome_name}_same_usearch.fa"
    family_file = f"{genome_name}_family_usearch.fa"
    
    process_fasta(df, trimmed_fasta, same_file, family_file, genome_name)
    
    print("Processing complete.")

if __name__ == "__main__":
    main()

