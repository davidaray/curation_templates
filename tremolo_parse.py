# Import necessary modules
import pandas as pd
import os
import argparse

# Use argparse to get input from the user
def GET_ARGS():
    parser = argparse.ArgumentParser(description="Process GC content and sequence length thresholds.")
    parser.add_argument("-i", "--input_list", type=str, required=True, help="List of species ids to find in the specified directory.")
    parser.add_argument("-d", "--directory", type=str, required=True, help="Directory to search for files.")
    parser.add_argument("-ocounts", "--output_counts", type=str, required=True, help="Ouput table file name.")
    
    # Get the arguments
    ARGS = parser.parse_args()
    INPUT_LIST = ARGS.input_list
    DIRECTORY = ARGS.directory
    OUTPUT_COUNTS = ARGS.output_counts
#    OUTPUT_RANKS = ARGS.output_ranks
    
    return INPUT_LIST, DIRECTORY, OUTPUT_COUNTS
    
def PROCESS_COUNTS(COUNTSFRAME, ID, DIRECTORY):
    PATH = os.path.join(DIRECTORY, ID + "_polymorphism_counts.tsv")
    if not os.path.exists(PATH):
        print(f"Warning: File {PATH} not found. Skipping {ID}.")
        return COUNTSFRAME  # Return without modifying the DataFrame
    INPUTCOUNTS = pd.read_csv(PATH, sep="\t", header=0)
    COUNTSFRAME[ID] = INPUTCOUNTS.iloc[:, 1]
    
if __name__ == "__main__": 

    # Get user inputs
    INPUT_LIST, DIRECTORY, OUTPUT_COUNTS = GET_ARGS()
    
    # Create an initial output file with the the first column
    DATA = {"TE_Class": ["Total", "DNA", "LINE", "LTR", "RC", "SINE"]}
    COUNTSFRAME = pd.DataFrame(DATA)

    # Parse input list and append counts to COUNTSFRAME
    with open(INPUT_LIST, "r") as LIST:
        for LINE in LIST:
            ID = LINE.strip()
            PROCESS_COUNTS(COUNTSFRAME, ID, DIRECTORY)
        
    COUNTSFRAME.to_csv(OUTPUT_COUNTS, sep='\t', index=False)  # index=False prevents writing row numbers

