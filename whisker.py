import pandas as pd
import glob
import gzip
import os
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

def parse_repeatmasker_file(file):
    """Parse a RepeatMasker .out.gz file and extract relevant information."""
    with gzip.open(file, 'rt') as f:
        # Read the file into a DataFrame, skipping the initial header lines.
        df = pd.read_csv(f, sep=r'\s+', header=None, skiprows=3,
                         names=['score', 'div.', 'del.', 'ins.', 'sequence', 'begin', 'end', 'left', 
                                'strand', 'matching', 'class/family', 'begin_R', 'end_R', 'left_R', 'ID'])
    return df

def parse_summary_file(file):
    """Extract total length from a RepeatMasker summary file."""
    with gzip.open(file, 'rt') as f:
        for line in f:
            if line.startswith("Total Length:"):
                total_length = int(line.split(":")[1].strip().split()[0])  # Get the length in bp
                return total_length
    return None

def load_mapping_file(mapping_file):
    """Load species to taxonomic family mapping from a file."""
    mapping_df = pd.read_csv(mapping_file, sep='\t')
    return dict(zip(mapping_df['Species_ID'], mapping_df['Taxonomic_Family']))

def process_repeatmasker_files(directory, species_to_family, total_length, output_prefix):
    """Process all RepeatMasker .out.gz files in a given directory and generate plots for each species."""
    all_files = glob.glob(os.path.join(directory, "*.out.gz"))

    for file in all_files:
        print(f"Processing {file}...")
        df = parse_repeatmasker_file(file)

        # Extract species ID from filename and map to taxonomic family.
        species_id = os.path.basename(file).split('.')[0].lstrip('m')[:6]
        species_family = species_to_family.get(species_id)

        if species_family:
            df['taxonomic_family'] = species_family

            # Bin divergence values.
            bins = [0, 5, 10, 15, 20, 25]
            df['div_bin'] = pd.cut(df['div.'], bins=bins)

            # Calculate counts and proportions for this species.
            count_df, proportion_df = calculate_counts_and_proportions(df, total_length)

            # Plot results for this species.
            plot_boxplots(count_df, proportion_df, f"{output_prefix}_{species_id}")

def calculate_counts_and_proportions(combined_df, total_length):
    """Calculate counts of insertions and proportions for each taxonomic family."""
    # Count insertions for each taxonomic family and divergence bin.
    count_df = combined_df.groupby(['taxonomic_family', 'div_bin']).size().reset_index(name='count')

    # Calculate proportion of bases occupied by each TE class.
    proportion_df = combined_df.groupby(['taxonomic_family', 'class/family']).agg(
        proportion=('score', 'sum')
    ).reset_index()

    # Normalize proportions using total_length
    proportion_df['proportion'] = proportion_df['proportion'] / total_length

    return count_df, proportion_df
    
def calculate_counts_and_proportions(combined_df, total_length):
    """Calculate counts of insertions and proportions for each taxonomic family."""
    # Count insertions for each taxonomic family and divergence bin.
    count_df = combined_df.groupby(['taxonomic_family', 'div_bin']).size().reset_index(name='count')

    # Calculate proportion of bases occupied by each TE class.
    proportion_df = combined_df.groupby(['taxonomic_family', 'class/family']).agg(
        proportion=('score', 'sum')
    ).reset_index()

    # Normalize proportions using total_length
    proportion_df['proportion'] = proportion_df['proportion'] / total_length

    return count_df, proportion_df

def plot_boxplots(count_df, proportion_df, output_prefix):
    """Create boxplots for counts and proportions."""
    # Plot 1: Boxplot of counts of insertions by taxonomic family.
    plt.figure(figsize=(12, 6))
    sns.boxplot(x='taxonomic_family', y='count', data=count_df)
    plt.xticks(rotation=45)
    plt.title('Box and Whisker Plot of TE Insertions by Taxonomic Family')
    plt.ylabel('Count of Insertions')
    plt.savefig(f"{output_prefix}_insertion_counts.png")
    plt.close()

    # Plot 2: Boxplot of proportions by taxonomic family.
    plt.figure(figsize=(12, 6))
    sns.boxplot(x='taxonomic_family', y='proportion', data=proportion_df)
    plt.xticks(rotation=45)
    plt.title('Box and Whisker Plot of TE Occupied Proportions by Taxonomic Family')
    plt.ylabel('Proportion of Bases Occupied')
    plt.savefig(f"{output_prefix}_proportions.png")
    plt.close()

def main():
    # Set up argument parsing.
    parser = argparse.ArgumentParser(description="Generate box and whisker plots for TE data.")
    parser.add_argument("-m", "--mapping", required=True, help="Mapping file of species to taxonomic families.")
    parser.add_argument("-d", "--directory", required=True, help="Directory containing RepeatMasker .out.gz files.")
    parser.add_argument("-o", "--output", required=True, help="Output prefix for the plots.")
    args = parser.parse_args()

    # Load species to family mapping.
    species_to_family = load_mapping_file(args.mapping)

    # Process RepeatMasker files and extract total length from the corresponding .summary.gz files.
    total_length = None
    summary_files = glob.glob(os.path.join(args.directory, "*.summary.gz"))

    if summary_files:
        total_length = parse_summary_file(summary_files[0])  # Assuming the same length for all

    if total_length is None:
        print("Total length not found in summary files.")
        return

    # Pass args.output to the processing function
    process_repeatmasker_files(args.directory, species_to_family, total_length, args.output)

if __name__ == "__main__":
    main()


"""
Starting over. Let's do this:

Starting with compressed out.gz files, I'd like to generate box and whisker plots for each family of bats. 

Input via argparse:
A mapping file of each species and it's corresponding taxonomic family will be provided, -m. 
A directory containing all .out.gz files will be provided, -d.
Calculate proportions of bases occupied using .summary.gz files provided in the same directory as the out.gz files.
Use this line from the .summary.gz files - Total Length: 1860258361 bp

For each taxonomic family, calculate box and whisker plots for the top five divergence bins (0-5, 5-10, 10-15, 15-20, 20-25). 

Two plots should be generated that encompass all Taxonomic families. 

Plot 1 - Y-axis should be the counts of insertions from the RepeatMasker .out.gz file. X- will be taxonomic families. Data should be aggregates for each taxonomic family.

Plot 2 - Y-axis should be number of proportion of bases occupied by each class of TE insertion (Classes to plot = DNA, RC, LINE, SINE, LTR, Unknown). X- will be taxonomic families. Data should be aggregates for each taxonomic family.

Example mapping data:
Taxonomic_Family	Binomial_Species_Name	Species_ID	Assembly_ID
Myzopodidae	Myzopoda aurita	MyzAur	mMyzAur1.1.pri
Emballonuridae	Rhynchonycteris naso	RhyNas	mRhyNas1.2.hap1
Emballonuridae	Saccopteryx bilineata	SacBil	mSacBil1.1.pri
Emballonuridae	Saccopteryx leptura	SacLep	mSacLep1.1.pri
Emballonuridae	Taphozous melanopogon	TapMel	mTapMel1.1.pri
Nycteridae	Nycteris thebaica	NycThe	mNycThe1.1.hap1

Use batches to conserve memory.

"""