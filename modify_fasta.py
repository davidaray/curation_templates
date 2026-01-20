import argparse

def parse_fasta(FILE):
    """Parse a fasta file and return a dictionary of headers and sequences."""
    SEQUENCES = {}
    HEADER = None
    with open(FILE, 'r') as F:
        for LINE in F:
            LINE = LINE.strip()
            if LINE.startswith(">"):
                HEADER = LINE
                SEQUENCES[HEADER] = []
            else:
                SEQUENCES[HEADER].append(LINE)
    for KEY in SEQUENCES:
        SEQUENCES[KEY] = ''.join(SEQUENCES[KEY])
    return SEQUENCES

def modify_headers(SEQUENCES, HEADERS_LIST, NEW_TEXT):
    """Modify headers by replacing the part after '#' with NEW_TEXT."""
    MODIFIED_SEQUENCES = {}
    for HEADER, SEQ in SEQUENCES.items():
        for OLD_HEADER in HEADERS_LIST:
            if OLD_HEADER in HEADER:
                NEW_HEADER = HEADER.split('#')[0] + f"#{NEW_TEXT}"
                MODIFIED_SEQUENCES[NEW_HEADER] = SEQ
                break
        else:
            MODIFIED_SEQUENCES[HEADER] = SEQ  # No change if HEADER not found
    return MODIFIED_SEQUENCES

def write_fasta(SEQUENCES, OUTPUT_FILE):
    """Write modified sequences to a new fasta file."""
    with open(OUTPUT_FILE, 'w') as F:
        for HEADER, SEQ in SEQUENCES.items():
            F.write(f"{HEADER}\n")
            F.write(f"{SEQ}\n")

def main():
    # Setup argparse
    PARSER = argparse.ArgumentParser(description="Modify fasta file headers")
    PARSER.add_argument("-f", "--fasta", required=True, help="Input fasta file")
    PARSER.add_argument("-i", "--header_ids", required=True, help="File with headers to modify")
    PARSER.add_argument("-n", "--newtext", required=True, help="New text for header modification")
    PARSER.add_argument("-o", "--output", required=True, help="Output fasta file")

    ARGS = PARSER.parse_args()

    # Load input files
    FASTA_FILE = ARGS.fasta
    HEADER_IDS_FILE = ARGS.header_ids
    NEW_TEXT = ARGS.newtext
    OUTPUT_FILE = ARGS.output

    # Read headers to modify
    with open(HEADER_IDS_FILE, 'r') as F:
        HEADERS_LIST = [LINE.strip() for LINE in F]

    # Parse fasta file
    SEQUENCES = parse_fasta(FASTA_FILE)

    # Modify headers
    MODIFIED_SEQUENCES = modify_headers(SEQUENCES, HEADERS_LIST, NEW_TEXT)

    # Save with user-defined output file name
    write_fasta(MODIFIED_SEQUENCES, OUTPUT_FILE)
    print(f"Modified fasta file saved as: {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
