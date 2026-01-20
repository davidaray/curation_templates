import os
import argparse
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description="Summarize TE content across species.")
    parser.add_argument("-m", "--mapping_file", required=True, help="Path to the mapping file (species_mapping_full.tsv).")
    parser.add_argument("-d", "--directory", required=True, help="Top-level directory containing the pie_data.tsv files.")
    parser.add_argument("-o", "--output", required=True, help="Output file for the summarized table.")
    return parser.parse_args()

def load_mapping(mapping_file):
    return pd.read_csv(mapping_file, sep="\t")

def extract_te_content(directory):
    te_data = {}
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith("pie_data.tsv"):
                species_id = os.path.basename(root).split("_")[0]
                file_path = os.path.join(root, file)
                data = pd.read_csv(file_path, sep="\t", index_col="Class")
                proportions = data["Proportion"].to_dict()
                te_data[species_id] = proportions
    return te_data

def summarize_data(mapping, te_data):
    columns = [
        "Taxonomic_Family", "Binomial_Species_Name", "Species_ID",
        "LINE_Proportion", "SINE_Proportion", "LTR_Proportion", 
        "DNA_Proportion", "RC_Proportion", "DIRS_Proportion", 
        "Unmasked_Proportion"
    ]
    summary = []
    for _, row in mapping.iterrows():
        species_id = row["Species_ID"]
        proportions = te_data.get(species_id, {})
        summary.append([
            row["Taxonomic_Family"], row["Binomial_Species_Name"], species_id,
            proportions.get("LINE", 0), proportions.get("SINE", 0),
            proportions.get("LTR", 0), proportions.get("DNA", 0),
            proportions.get("RC", 0), proportions.get("DIRS", 0),
            proportions.get("Unmasked", 0)
        ])
    return pd.DataFrame(summary, columns=columns)

def main():
    args = parse_arguments()
    
    # Load the mapping file
    mapping = load_mapping(args.mapping_file)
    
    # Extract TE content from all pie_data.tsv files
    te_data = extract_te_content(args.directory)
    
    # Summarize the data
    summary_table = summarize_data(mapping, te_data)
    
    # Save the output
    summary_table.to_csv(args.output, sep="\t", index=False)
    print(f"Summarized table saved to {args.output}")

if __name__ == "__main__":
    main()

