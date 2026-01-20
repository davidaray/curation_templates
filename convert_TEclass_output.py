#!/usr/bin/env python3
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Process TE classification table")
    parser.add_argument("-i", "--input", required=True, help="Input TSV file")
    parser.add_argument("-o1", "--output1", required=True, help="Output TSV (raw formatted)")
    parser.add_argument("-o2", "--output2", required=True, help="Output TSV (mapped + filtered)")
    args = parser.parse_args()

    # Load table
    df = pd.read_csv(args.input, sep="\t")

    # Keep only needed columns
    df = df[["name", "order", "class", "probability"]].copy()

    # Remove anything after whitespace in 'name'
    df["name"] = df["name"].str.split().str[0]

    # Save raw formatted table
    df.to_csv(args.output1, sep="\t", index=False, header=["Name", "Order", "Class", "Probability"])

    # Mapping dictionary
    mapping = {
        ("Class_I-LINE", "L1_L2"): "LINE/L1",
        ("Class_I-LTR", "ERV"): "LTR/ERV",
        ("Class_I", "SINE"): "SINE/tRNA",
        ("Unknown", "Unknown"): "Unknown/Unknown",
        ("Class_I-LINE", "RTE"): "LINE/RTE",
        ("Class_II-TIR", "TcMar"): "DNA/TcMariner",
        ("Class_I-LTR", "Copia"): "LTR/Copia",
        ("Class_I-LTR", "Gypsy"): "LTR/Gypsy",
        ("Class_II-TIR", "hAT"): "DNA/hAT",
        ("Class_I-LINE", "Jockey"): "LINE/Jockey"
    }

    # Replace Order/Class with mapped values
    df["SuperClass"] = df.apply(lambda x: mapping.get((x["order"], x["class"]), f"{x['order']}/{x['class']}"), axis=1)

    # Reorder columns
    df = df[["name", "SuperClass", "probability"]]

    # Filter rows with probability >= 0.7
    df_filtered = df[df["probability"] >= 0.7].copy()

    # Save mapped + filtered table
    df_filtered.to_csv(args.output2, sep="\t", index=False, header=["Name", "SuperClass", "Probability"])

if __name__ == "__main__":
    main()

