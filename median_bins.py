import pandas as pd
import argparse
import os

def find_median_bin(data, column):
    """
    Calculate the median bin location for a given column in the dataset.
    :param data: DataFrame containing the data.
    :param column: Name of the column for which to calculate the median bin.
    :return: Median bin index.
    """
    total = data[column].sum()
    cumulative = data[column].cumsum()
    median_bin = data[cumulative >= total / 2].iloc[0]["Bin"]
    return median_bin

def extract_label_from_filename(filename):
    """
    Extract the label from the input filename by removing the file extension and trailing identifiers.
    :param filename: Input filename.
    :return: Extracted label.
    """
    basename = os.path.basename(filename)
    label = basename.split("_stacked_bar_data")[0]
    return label

def calculate_median_bins(input_file, output_file):
    """
    Calculate median bins for all TE classes and save results to an output file in transposed format.
    :param input_file: Path to the input data file.
    :param output_file: Path to the output results file.
    """
    # Read input data
    df = pd.read_csv(input_file, delim_whitespace=True)

    # Ensure "Bin" is numeric
    df["Bin"] = pd.to_numeric(df["Bin"], errors="coerce")

    # Calculate median bins for each class
    median_bins = {col: find_median_bin(df, col) for col in df.columns if col != "Bin"}

    # Extract label from the input filename
    label = extract_label_from_filename(input_file)

    # Convert the median bins to a DataFrame
    result_df = pd.DataFrame(median_bins, index=[label])

    # Write the transposed DataFrame to the output file
    result_df.to_csv(output_file, sep="\t", header=True, index=True)
    print(f"Results saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Calculate median bin locations for TE classes.")
    parser.add_argument("-i", "--input", required=True, help="Path to input data file")
    parser.add_argument("-o", "--output", required=True, help="Path to output results file")
    args = parser.parse_args()

    calculate_median_bins(args.input, args.output)

if __name__ == "__main__":
    main()
