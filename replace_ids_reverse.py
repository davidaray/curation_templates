import argparse
import pandas as pd
import re

def load_mapping(mapping_file):
    """Load the mapping file into a dictionary with Binomial_Species_Name as the key and Old_ID as the value."""
    mapping = {}
    with open(mapping_file, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    old_id = parts[4]  # Old_ID column
                    species_name = parts[1] # Binomial_Species_Name column
                    mapping[species_name] = old_id  # Reverse mapping
    return mapping

def replace_names_in_file(input_file, mapping):
    """Replace Binomial_Species_Name occurrences in the file with Old_ID."""
    with open(input_file, 'r') as f:
        content = f.read()
    
    for species_name, old_id in mapping.items():
        content = re.sub(rf'\b{species_name}\b', old_id, content)
    
    with open(input_file, 'w') as f:
        f.write(content)

def main():
    parser = argparse.ArgumentParser(description="Replace Binomial_Species_Names in a file with corresponding Old_IDs.")
    parser.add_argument("-m", "--mapping", required=True, help="Path to the mapping file")
    parser.add_argument("-f", "--file", required=True, help="Path to the file to be modified")
    args = parser.parse_args()
    
    mapping = load_mapping(args.mapping)
    replace_names_in_file(args.file, mapping)
    print("Reverse replacement complete.")

if __name__ == "__main__":
    main()
