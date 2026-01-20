import os
import argparse
import pandas as pd

def process_files(input_folder, output_file, filter_value):
    # Collect all files in the input folder
    all_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if os.path.isfile(os.path.join(input_folder, f))]
    
    # Initialize an empty dataframe to store concatenated data
    concatenated_data = pd.DataFrame()
    
    for file in all_files:
        # Read each file into a pandas dataframe
        df = pd.read_csv(file, delim_whitespace=True, header=None, skiprows=1)
        concatenated_data = pd.concat([concatenated_data, df], ignore_index=True)
    
    # Filter rows where the value in column 5 (index 4) is >= filter_value
    filtered_data = concatenated_data[concatenated_data[4] >= filter_value]
    
    # Select the first four columns
    result = filtered_data[[0, 1, 2, 3]]
    
    # Save to output file
    result.to_csv(output_file, sep='\t', header=False, index=False)

def main():
    # Set up argparse to handle command line arguments
    parser = argparse.ArgumentParser(description="Process files, filter, and save output.")
    parser.add_argument('-if', '--input_folder', required=True, help="Input folder containing the files")
    parser.add_argument('-o', '--output_file', required=True, help="Output file name")
    parser.add_argument('-f', '--filter_value', required=True, type=float, help="Filter value for column 5 (e.g., 0.7)")
    
    args = parser.parse_args()
    
    # Process the files and create the output file
    process_files(args.input_folder, args.output_file, args.filter_value)

if __name__ == '__main__':
    main()
