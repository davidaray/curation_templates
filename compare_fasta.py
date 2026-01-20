import argparse
from Bio import SeqIO

def compare_fasta_files(file1, file2):
    # Read the sequences from the first FASTA file into a dictionary
    sequences1 = SeqIO.to_dict(SeqIO.parse(file1, "fasta"))
    
    # Read the sequences from the second FASTA file into a dictionary
    sequences2 = SeqIO.to_dict(SeqIO.parse(file2, "fasta"))
    
    # Find sequences that are in file1 but not in file2
    unique_to_file1 = set(sequences1.keys()) - set(sequences2.keys())
    
    # Find sequences that are in file2 but not in file1
    unique_to_file2 = set(sequences2.keys()) - set(sequences1.keys())
    
    # Find sequences that are in both files but differ in sequence identity
    differing_sequences = []
    for key in set(sequences1.keys()).intersection(set(sequences2.keys())):
        if str(sequences1[key].seq) != str(sequences2[key].seq):
            differing_sequences.append(key)
    
    return unique_to_file1, unique_to_file2, differing_sequences

def main():
    parser = argparse.ArgumentParser(description="Identify differences in sequence identity between two FASTA files.")
    parser.add_argument('-f1', '--fasta1', required=True, help="Path to the first FASTA file.")
    parser.add_argument('-f2', '--fasta2', required=True, help="Path to the second FASTA file.")
    
    args = parser.parse_args()
    
    unique_to_file1, unique_to_file2, differing_sequences = compare_fasta_files(args.fasta1, args.fasta2)
    
    if unique_to_file1:
        print(f"Sequences unique to {args.fasta1}:")
        for seq_id in unique_to_file1:
            print(seq_id)
    
    if unique_to_file2:
        print(f"Sequences unique to {args.fasta2}:")
        for seq_id in unique_to_file2:
            print(seq_id)
    
    if differing_sequences:
        print("Sequences present in both files but with different identities:")
        for seq_id in differing_sequences:
            print(seq_id)

if __name__ == "__main__":
    main()

