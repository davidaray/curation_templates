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
        
        if classification == 'same':
            # Write the sequence to the "known elements" file and remove it from new_seqs
            SeqIO.write(new_seqs[fasta_header_new], known_elements_output, 'fasta')
            del new_seqs[fasta_header_new]
        else:
            # Write the sequence to the concatenated output file
            SeqIO.write(new_seqs[fasta_header_new], concatenated_output, 'fasta')
    
    # Write remaining sequences to the concatenated output file
    for seq in new_seqs.values():
        SeqIO.write(seq, concatenated_output, 'fasta')
    for seq in known_seqs.values():
        SeqIO.write(seq, concatenated_output, 'fasta')

    known_elements_output.close()
    concatenated_output.close()

