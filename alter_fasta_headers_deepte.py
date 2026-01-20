import argparse

def load_mapping(MAPPING_FILE):
    """Loads the mapping of column A to column B from the input file."""
    MAPPING_DICT = {}
    with open(MAPPING_FILE, 'r') as F:
        for LINE in F:
            COL_A, COL_B = LINE.strip().split('\t')
            MAPPING_DICT[COL_A] = COL_B
    return MAPPING_DICT

def modify_fasta_headers(INPUT_FASTA, MAPPING_DICT, OUTPUT_FASTA):
    """Modifies the headers of the input FASTA file based on the mapping and writes to the output FASTA file."""
    with open(INPUT_FASTA, 'r') as INFILE, open(OUTPUT_FASTA, 'w') as OUTFILE:
        for LINE in INFILE:
            if LINE.startswith('>'):
                # Split the header to isolate the part after '#'
                HEADER = LINE.strip()
                PARTS = HEADER.split('#')
                if len(PARTS) > 1:
                    # Get the part after '#' and extract class information
                    CLASS_INFO = PARTS[1].split('__')[-1]
                    # Replace the part after '#' with the mapping value if found
                    if CLASS_INFO in MAPPING_DICT:
                        NEW_HEADER = PARTS[0] + '#' + MAPPING_DICT[CLASS_INFO]
                        OUTFILE.write(NEW_HEADER + '\n')
                    else:
                        # If no mapping is found, keep the original header
                        OUTFILE.write(HEADER + '\n')
                else:
                    # In case no '#' is present, keep the original header
                    OUTFILE.write(HEADER + '\n')
            else:
                # Write the sequence lines unchanged
                OUTFILE.write(LINE)

def main():
    PARSER = argparse.ArgumentParser(description="Modify FASTA headers based on a mapping file.")
    PARSER.add_argument("-f", "--fasta", required=True, help="Input FASTA file")
    PARSER.add_argument("-m", "--mapping", required=True, help="Input mapping file with two columns")
    
    ARGS = PARSER.parse_args()

    INPUT_FASTA = ARGS.fasta
    MAPPING_FILE = ARGS.mapping
    OUTPUT_FASTA = INPUT_FASTA.rsplit('.', 1)[0] + '_deep_te.fa'

    # Load the mapping dictionary from the input file
    MAPPING_DICT = load_mapping(MAPPING_FILE)

    # Modify the FASTA headers and write to a new file
    modify_fasta_headers(INPUT_FASTA, MAPPING_DICT, OUTPUT_FASTA)
    print(f"Modified FASTA saved to: {OUTPUT_FASTA}")

if __name__ == "__main__":
    main()

