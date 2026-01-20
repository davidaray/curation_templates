##Create a python script, hite_usearch_final.py, that will take in 3 arguments using argparse.
#1. a "known library" file, -m
#2. a "new library" file called ${ID}_final_library.fa, -n
#3. a .tsv file, -t

##The .tsv file should have five columns, A-E.
#A=fasta header from "new library"
#B=fasta header from "known library"
#C=percent identity
#D=length of sequence from new library
#E=length of sequence from known library

##The python script should then:
##Add a new column, F, to the .tsv file. The value in column F should be the value in D divided by E.
##Add a new column, G, to the .tsv file. 
#	If column C>95 and D>80 and F is between 0.8 and 1.2, enter 'same' in column G.
#	If column C>80 but less than 95, and D>80, and F is between 0.8 and 1.2, enter 'family' in column G. 
#	If column C>60 but less than 80, and D>80, and F is between 0.8 and 1.2, enter 'class' in column G. 
##Save the file with the new columns as ${ID}_usearch_final.tsv.

##Process ${ID}_usearch_final.tsv as follows:
#For any rows where column G is 'same', remove the associated sequence from ${ID}_final_library.fa and save #to a file called "${ID)_known_elements_from_final_search.fa". For all other sequences, add them to the #"known library" file and save the new file as "../concatenated_library_${ID}.fa".

import argparse
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

def process_tsv(tsv_file):
    # Read the TSV file into a DataFrame
    df = pd.read_csv(tsv_file, sep='\t', header=None)
    df.columns = ['A', 'B', 'C', 'D', 'E']
    
    # Add column F: D divided by E
    df['F'] = df['D'] / df['E']
    
    # Add column G with the conditions specified
    conditions = [
        (df['C'] > 95) & (df['D'] > 80) & (df['F'].between(0.8, 1.2)),
        (df['C'] > 80) & (df['C'] <= 95) & (df['D'] > 80) & (df['F'].between(0.8, 1.2)),
        (df['C'] > 60) & (df['C'] <= 80) & (df['D'] > 80) & (df['F'].between(0.8, 1.2))
    ]
    choices = ['same', 'family', 'class']
    df['G'] = pd.Series(pd.cut(df['C'], bins=[0, 60, 80, 95, 100], labels=['', 'class', 'family', 'same']))[
        df['F'].between(0.8, 1.2) & (df['D'] > 80)
    ].fillna('')[df['C'] > 60]

    # Save the updated TSV file
    new_tsv_file = tsv_file.replace('.tsv', '_usearch_final.tsv')
    df.to_csv(new_tsv_file, sep='\t', index=False)

    return df, new_tsv_file

def process_fasta_files(new_library, known_library, tsv_data, output_id):
    # Read the sequences from the new library
    new_seqs = SeqIO.to_dict(SeqIO.parse(new_library, 'fasta'))
    
    # Handle duplicate keys in new_seqs
    new_seq_keys = defaultdict(int)
    for key in list(new_seqs.keys()):
        new_seq_keys[key] += 1
        if new_seq_keys[key] > 1:
            new_key = f"{key}_{new_seq_keys[key]}"
            new_seqs[new_key] = new_seqs.pop(key)
    
    # Read the sequences from the known library and handle duplicates
    known_seq_keys = defaultdict(int)
    known_seqs = {}
    for record in SeqIO.parse(known_library, 'fasta'):
        key = record.id
        known_seq_keys[key] += 1
        if known_seq_keys[key] > 1:
            key = f"{key}_{known_seq_keys[key]}"
        known_seqs[key] = record

    # Open output files
    known_elements_output = open(f"{output_id}_known_elements_from_final_search.fa", 'w')
    concatenated_output = open(f"../concatenated_library_{output_id}.fa", 'w')
    
    # Process each row in the TSV file
    for index, row in tsv_data.iterrows():
        fasta_header_new = row['A']
        fasta_header_known = row['B']
        classification = row['G']
        
        if fasta_header_new not in new_seqs:
            print(f"Warning: {fasta_header_new} not found in new_seqs. It was moved to the duplicates file. Skipping this entry.")
            continue
        
        if classification == 'same':
            # Write the sequence to the "known elements" file and remove it from new_seqs
            SeqIO.write(new_seqs[fasta_header_new], known_elements_output, 'fasta')
            del new_seqs[fasta_header_new]
        else:
            # Write the sequence to the concatenated output file
            SeqIO.write(new_seqs[fasta_header_new], concatenated_output, 'fasta')
    
    # Write remaining sequences to the concatenated output file
    for seq in known_seqs.values():
        SeqIO.write(seq, concatenated_output, 'fasta')

    known_elements_output.close()
    concatenated_output.close()
    concatenated_output_file = f"../concatenated_library_{output_id}.fa"
    print('Output file is ' + concatenated_output_file)

def main():
    parser = argparse.ArgumentParser(description="Process usearch final results.")
    parser.add_argument('-m', '--known_library', required=True, help="Path to the known library fasta file.")
    parser.add_argument('-n', '--new_library', required=True, help="Path to the new library fasta file.")
    parser.add_argument('-t', '--tsv_file', required=True, help="Path to the TSV file.")
    
    args = parser.parse_args()

    # Extract the ID from the new library filename
    output_id = args.new_library.split('_')[0]

    # Process the TSV file and save the new version with additional columns
    tsv_data, new_tsv_file = process_tsv(args.tsv_file)

    # Process the fasta files and create output files
    process_fasta_files(args.new_library, args.known_library, tsv_data, output_id)

if __name__ == "__main__":
    main()
