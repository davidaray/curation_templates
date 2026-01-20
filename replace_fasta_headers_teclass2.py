import argparse
import re

def load_replacement_file(replacement_file):
    """Loads the replacement data from the input file and creates a dictionary for fast lookup."""
    replacement_dict = {}
    with open(replacement_file, 'r') as f:
        for line in f:
            cols = line.strip().split()
            if len(cols) >= 4:
                header_part = cols[0].strip()
                col3 = cols[2].strip()
                col4 = cols[3].strip()
                replacement_dict[header_part] = f"{col3}/{col4}"
            else:
                print(f"Skipping line in replacement file (not enough columns): {line.strip()}")
    return replacement_dict

def apply_final_fixes(header):
    """Applies final fixes to the header."""
    # Remove "Class_I-" from the header
    header = header.replace("Class_I-", "")
    
    # Replace "Class_I/SINE" with "SINE/tRNA"
    header = header.replace("Class_I/SINE", "SINE/tRNA")
    
    # Remove "Class_II-" from the header
    header = header.replace("Class_II-", "")
    
    # Replace "TIR/" with "DNA/"
    header = header.replace("TIR/", "DNA/")
    
    #Replace "TcMar" with "TcMariner", only in cases where TcMar is entire word
    header = re.sub(r'\bTcMar\b', 'TcMariner', header)
        
    
    return header

def modify_fasta_headers(fasta_file, replacement_dict, output_file):
    """Modifies the headers of the input FASTA file based on the replacement file and applies final fixes."""
    with open(fasta_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                header = line.strip()
                # Remove '>' from header
                header_no_gt = header.lstrip('>')
                # Split on '#' if present
                if '#' in header_no_gt:
                    header_before_hash, _ = header_no_gt.split('#', 1)
                else:
                    header_before_hash = header_no_gt
                # Check if the part before '#' exists in the replacement dictionary
                if header_before_hash in replacement_dict:
                    new_classification = replacement_dict[header_before_hash]
                    new_header = f">{header_before_hash}#{new_classification}"
                    # Apply final fixes to the new header
                    new_header = apply_final_fixes(new_header)
                    print(f"Header found: {header}")
                    print(f"New header: {new_header}")
                    outfile.write(new_header + '\n')
                else:
                    print(f"Header not found: {header}")
                    outfile.write(line)
            else:
                # Write the sequence lines unchanged
                outfile.write(line)

def main():
    parser = argparse.ArgumentParser(description="Modify FASTA headers based on a replacement file.")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file with DNA sequences")
    parser.add_argument("-r", "--replacement", required=True, help="Replacement file with headers and new classifications")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    
    args = parser.parse_args()

    # Load the replacement data from the replacement file
    replacement_dict = load_replacement_file(args.replacement)

    # Modify the FASTA headers based on the replacement file
    modify_fasta_headers(args.fasta, replacement_dict, args.output)

    print(f"Modified FASTA saved to: {args.output}")

if __name__ == "__main__":
    main()
