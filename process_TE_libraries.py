import argparse
import os
import subprocess
import shutil

# Hardcoded paths for MMseqs2
MMSEQS_PATH = '/lustre/work/daray/software/MMseqs2/build/bin/mmseqs'
TMP_DIR = '/lustre/scratch/daray/bat1k_TE_analyses/mmseqs_tmp'

# Number of threads to use
NUM_THREADS = 36

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process TE libraries and identify duplicates using MMseqs2.')
    parser.add_argument('-i', '--id_list', required=True, help='Path to the text list of IDs.')
    parser.add_argument('-l', '--library', required=True, help='Path to the known TEs library (mammals.plus.covid_bats2.04072022.fa).')
    parser.add_argument('-o', '--output_dir', required=True, help='Directory to store the output files.')
    return parser.parse_args()

def concatenate_files(id_list, library, output_dir):
    concatenated_file = os.path.join(output_dir, 'concatenated.fa')
    with open(concatenated_file, 'w') as outfile:
        # Write the contents of the known TEs library to the concatenated file
        with open(library, 'r') as infile:
            shutil.copyfileobj(infile, outfile)
        
        # Iterate over the IDs and append the contents of each _final_library.fa
        for id_ in id_list:
            final_library_file = os.path.join(id_, f'{id_}_final_library.fa')
            if os.path.exists(final_library_file):
                with open(final_library_file, 'r') as infile:
                    shutil.copyfileobj(infile, outfile)
            else:
                print(f"Warning: {final_library_file} not found. Skipping.")

    return concatenated_file

def run_mmseqs(concatenated_file, output_dir):
    # Create necessary directories
    os.makedirs(TMP_DIR, exist_ok=True)
    db = os.path.join(output_dir, 'db')
    result = os.path.join(output_dir, 'result')
    tsv_file = os.path.join(output_dir, 'mmseqs_results.tsv')

    # Convert FASTA to MMseqs2 database
    subprocess.run([MMSEQS_PATH, 'createdb', concatenated_file, db], check=True)

    # Run MMseqs2 search with 36 threads
    subprocess.run([
        MMSEQS_PATH, 'search', db, db, result, TMP_DIR,
        '--min-seq-id', '0.95',  # 95% identity
        '--cov-mode', '0',       # Coverage mode
        '--max-seqs', '100',
        '--search-type', '3',
        '--threads', str(NUM_THREADS),
        '--local-tmp',  '/lustre/scratch/daray/bat1k_TE_analyses/mmseqs_tmp',
        '-v', '3',
        '-s', '7'
    ], check=True)

    # Convert MMseqs2 result to TSV
    subprocess.run([MMSEQS_PATH, 'convertalis', db, db, result, tsv_file], check=True)

    return tsv_file

def filter_sequences(tsv_file, unique_file, duplicates_file, library):
    # Read the sequences from the known TEs library
    known_tes = set()
    with open(library, 'r') as lib:
        seq_id = None
        for line in lib:
            if line.startswith('>'):
                seq_id = line.split()[0][1:]  # Extract sequence ID without '>'
                known_tes.add(seq_id)
    
    # Track seen sequences and lengths
    seen_sequences = {}
    with open(tsv_file, 'r') as tsv, open(unique_file, 'w') as unique_out, open(duplicates_file, 'w') as duplicates_out:
        for line in tsv:
            query, target, pident, alnlen, qlen, tlen, *_ = line.strip().split('\t')
            qlen, tlen = int(qlen), int(tlen)
            length_ratio = min(qlen / tlen, tlen / qlen)

            if 0.8 <= length_ratio <= 1.2:
                if query in seen_sequences or target in seen_sequences:
                    duplicates_out.write(f">{query}\n")
                else:
                    if target in known_tes:
                        unique_out.write(f">{target}\n")
                        seen_sequences[target] = tlen
                    else:
                        unique_out.write(f">{query}\n")
                        seen_sequences[query] = qlen

def main():
    args = parse_arguments()

    # Ensure the output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # Read the list of IDs from the input file
    with open(args.id_list, 'r') as f:
        id_list = [line.strip() for line in f if line.strip()]

    # Concatenate the _final_library.fa files with the known TEs library
    concatenated_file = concatenate_files(id_list, args.library, args.output_dir)

    # Run MMseqs2 to identify and cluster sequences based on similarity
    tsv_file = run_mmseqs(concatenated_file, args.output_dir)

    # Filter sequences, keeping those from the known TEs library and removing duplicates
    filter_sequences(tsv_file, os.path.join(args.output_dir, 'unique_sequences.fa'), os.path.join(args.output_dir, 'duplicates.fa'), args.library)

    print("Processing complete.")

if __name__ == '__main__':
    main()
