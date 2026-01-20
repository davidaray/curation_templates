import os
import argparse
import pandas as pd

def process_file(input_file, output_file, filter_value):
    """
    Process a single TSV file, filter by column 3 value, and save output.
    
    Args:
        input_file: Path to input TSV file
        output_file: Path to output TSV file
        filter_value: Minimum value for filtering column 3
    """
    # Read the file into a pandas dataframe
    df = pd.read_csv(input_file, delim_whitespace=True, header=None, skiprows=1)
    
    # Filter rows where the value in column 3 (index 2) is >= filter_value
    filtered_data = df[df[2] >= filter_value]
    
    # Select the first three columns
    result = filtered_data[[0, 1, 2]]
    
    # Save to output file
    result.to_csv(output_file, sep='\t', header=False, index=False)
    
    print(f"Processed {len(df)} rows")
    print(f"Filtered to {len(result)} rows where column 5 >= {filter_value}")
    print(f"Output saved to: {output_file}")

def main():
    # Set up argparse to handle command line arguments
    parser = argparse.ArgumentParser(description="Process TSV file, filter by column 5, and save output.")
    parser.add_argument('-i', '--input_file', required=True, help="Input TSV file")
    parser.add_argument('-o', '--output_file', required=True, help="Output file name")
    parser.add_argument('-f', '--filter_value', required=True, type=float, help="Filter value for column 3 (e.g., 0.7)")
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.isfile(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found.")
        return
    
    # Process the file and create the output file
    process_file(args.input_file, args.output_file, args.filter_value)

if __name__ == '__main__':
    main()
