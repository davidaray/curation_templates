import argparse
import os
import pandas as pd
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description='Process genome files.')
    parser.add_argument('-u', '--usearch', required=True, help='Path to ${GENOME}_usearch_60_hits.tsv')
    parser.add_argument('-c', '--cons', required=True, help='Path to ${GENOME}_$(TETYPE)_collection-_rep.fa')
    parser.add_argument('-m', '--mammals', required=True, help='Path to library of known elements')
    parser.add_argument('-t', '--tetype', required=True, help='TE type label for output files')
    return parser.parse_args()

def process_usearch_file(usearch_file):
    df = pd.read_csv(usearch_file, sep='\t', header=None)
    df.columns = ['A', 'B', 'C', 'D', 'E']
    df['F'] = df['D'] / df['E']
    df['G'] = df.apply(lambda row: 'same' if row['C'] > 95 and row['D'] > 80 and 0.8 <= row['F'] <= 1.2 
                       else ('family' if 80 <= row['C'] <= 95 and row['D'] > 80 and 0.8 <= row['F'] <= 1.2 
                             else ('class' if 60 <= row['C'] < 80 and row['D'] > 80 and 0.8 <= row['F'] <= 1.2 else '')), axis=1)
    df.to_csv(usearch_file, sep='\t', index=False, header=False)
    return df

def copy_and_trim_fasta(cons_file, trimmed_fasta):
    SeqIO.write(SeqIO.parse(cons_file, "fasta"), trimmed_fasta, "fasta")

def process_fasta(df, cons_file, same_file, family_file, class_file, genome_name, tetype):
    trimmed_seqs = {record.id: record for record in SeqIO.parse(cons_file, "fasta")}
    
    same_seqs = []
    family_seqs = []
    class_seqs = []
    
    for _, row in df.iterrows():
        if row['G'] == 'same':
            if row['A'] in trimmed_seqs:
                record = trimmed_seqs.pop(row['A'])
                same_seqs.append(record)
                print(f"EXACT MATCH: Replaced rediscovered {row['A']} with {record.id}")
        elif row['G'] == 'family':
            if row['A'] in trimmed_seqs:
                record = trimmed_seqs.pop(row['A'])
                new_header = row['B'].replace('#', f"_{genome_name}#")
                record.id = new_header
                record.description = ''
                family_seqs.append(record)
                print(f"FAMILY MATCH: Replaced header for {row['A']} with {record.id}")
        elif row['G'] == 'class':
            if row['A'] in trimmed_seqs:
                record = trimmed_seqs.pop(row['A'])
                new_header = row['B'].replace('#', f"_{genome_name}#")
                record.id = new_header
                record.description = ''
                class_seqs.append(record)
                print(f"CLASS MATCH: Replaced header for {row['A']} with {record.id}")
    
    SeqIO.write(same_seqs, same_file, "fasta")
    SeqIO.write(family_seqs, family_file, "fasta")
    SeqIO.write(class_seqs, class_file, "fasta")
    SeqIO.write(trimmed_seqs.values(), cons_file, "fasta")

def rename_and_classify_fasta(trimmed_fasta, renamed_fasta, genome_name):
    def classify_header(header):
        if "Homology_Non_LTR" in header or "Denovo_Non_LTR" in header:
            return "#NonLTR/Unknown"
        elif "Helitron" in header:
            return "#RC/Unknown"
        elif "LTR_" in header:
            if "_INT" in header:
                return "#LTR/Unknown"
            elif "_LTR" in header:
                return "#LTR/Unknown"
        elif "TIR" in header:
            return "#DNA/Unknown"
        else:
            return "#Unknown"

    renamed_seqs = []
    for record in SeqIO.parse(trimmed_fasta, "fasta"):
        classification = classify_header(record.id)
        new_id = f"{genome_name}_{record.id}{classification}"
        record.id = new_id
        record.description = ''
        renamed_seqs.append(record)
    
    SeqIO.write(renamed_seqs, renamed_fasta, "fasta")

def main():
    args = parse_args()
    
    genome_name = os.path.basename(args.usearch).split('_')[0]
    tetype = args.tetype
    
    df = process_usearch_file(args.usearch)
    
    trimmed_fasta = f"{genome_name}_{tetype}_collection.cons.trimmed.fa"
    copy_and_trim_fasta(args.cons, trimmed_fasta)
    
    same_file = f"{genome_name}_{tetype}_same_usearch.fa"
    family_file = f"{genome_name}_{tetype}_family_usearch.fa"
    class_file = f"{genome_name}_{tetype}_class_usearch.fa"
    
    process_fasta(df, trimmed_fasta, same_file, family_file, class_file, genome_name, tetype)
    
    renamed_fasta = f"{genome_name}_{tetype}_cons.renamed.fa"
    rename_and_classify_fasta(trimmed_fasta, renamed_fasta, genome_name)
    
    print("Processing complete.")

if __name__ == "__main__":
    main()
