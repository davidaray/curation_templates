import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

def load_data(file_path):
    """
    Load the input TSV file into a pandas DataFrame.
    """
    return pd.read_csv(file_path, sep="\t")

def plot_boxplots(data, output_file):
    """
    Create boxplots of TE class proportions by taxonomic family.
    """
    # Melt the DataFrame for easier plotting
    df_melted = data.melt(
        id_vars="Taxonomic_Family", 
        value_vars=["LINE_Proportion", "SINE_Proportion", "LTR_Proportion", "DNA_Proportion", "RC_Proportion"],
        var_name="TE_Class", 
        value_name="Proportion"
    )

    # Create boxplots
    plt.figure(figsize=(12, 8))
    sns.boxplot(data=df_melted, x="TE_Class", y="Proportion", hue="Taxonomic_Family")
    plt.xticks(rotation=45)
    plt.title("TE Class Proportions by Taxonomic Family", fontsize=16)
    plt.ylabel("Proportion")
    plt.xlabel("TE Class")
    plt.legend(title="Taxonomic Family", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    # Save the plot to the specified output file
    plt.savefig(output_file)
    plt.close()

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Generate boxplots of TE class proportions by taxonomic family.")
    parser.add_argument("-i", "--input", required=True, help="Input TSV file with TE class proportions.")
    parser.add_argument("-o", "--output", required=True, help="Output file for the boxplot image (e.g., 'output.png').")
    args = parser.parse_args()

    # Load the data and generate the plot
    data = load_data(args.input)
    plot_boxplots(data, args.output)

if __name__ == "__main__":
    main()

