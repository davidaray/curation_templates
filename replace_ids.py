import argparse
import pandas as pd
import re

def load_mapping(mapping_file):
    """Load the mapping file into a dictionary."""
    mapping = {}
    with open(mapping_file, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    old_id = parts[4]  # Old_ID column
                    species_name = parts[1].replace(" ", "_")  # Binomial_Species_Name column
                    mapping[old_id] = species_name
    return mapping

def replace_ids_in_file(input_file, mapping):
    """Replace Old_ID occurrences in the file with Binomial_Species_Name."""
    with open(input_file, 'r') as f:
        content = f.read()
    
    for old_id, species_name in mapping.items():
        content = re.sub(rf'\b{old_id}\b', species_name, content)
    
    with open(input_file, 'w') as f:
        f.write(content)

def main():
    parser = argparse.ArgumentParser(description="Replace Old_IDs in a file with corresponding Binomial_Species_Names.")
    parser.add_argument("-m", "--mapping", required=True, help="Path to the mapping file")
    parser.add_argument("-f", "--file", required=True, help="Path to the file to be modified")
    args = parser.parse_args()
    
    mapping = load_mapping(args.mapping)
    replace_ids_in_file(args.file, mapping)
    print("Replacement complete.")

if __name__ == "__main__":
    main()
