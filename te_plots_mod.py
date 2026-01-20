import argparse
import pandas as pd
# Use non-interactive Agg backend
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gzip
import re
import multiprocessing as mp

# Define color dictionary for TE classes
class_colors = {
    'DNA': '#1f77b4', 'DIRS': '#ff7f0e', 'LINE': '#2ca02c', 'LTR': '#d62728', 
    'RC': '#9467bd', 'SINE': '#8c564b', 'NonLTR': '#e377c2', 'Satellite': '#7f7f7f', 
    'Unknown': '#bcbd22', 'Other': '#17becf', 'Non-TE': '#000000'
}

# Function to parse the .out file (supporting .gz)
def parse_repeatmasker_outfile(out_file):
    columns = ['SW_score', 'perc_div', 'perc_del', 'perc_ins', 'query_sequence', 'begin', 'end', 'left', 
               'strand', 'matching_repeat', 'repeat_class', 'repeat_pos_begin', 'repeat_pos_end', 'repeat_left', 'ID']
    
    if out_file.endswith('.gz'):
        with gzip.open(out_file, 'rt') as f:
            df = pd.read_csv(f, delim_whitespace=True, skiprows=3, names=columns, engine='python')
    else:
        df = pd.read_csv(out_file, delim_whitespace=True, skiprows=3, names=columns, engine='python')

    df['TE_class'] = df['repeat_class'].apply(lambda x: x.split('/')[0] if isinstance(x, str) and '/' in x else x)
    df['TE_class'] = df['TE_class'].apply(lambda x: 'Satellite' if x == 'Simple_repeat' else x)
    df['TE_class'] = df['TE_class'].apply(lambda x: 'Unknown' if not isinstance(x, str) else x)
    df['insertion_length'] = df['end'] - df['begin']
    
    return df

# Function to parse the .bed file
def parse_bed_file(bed_file):
    columns = ['scaffold', 'start', 'end', 'repeat_name', 'length', 'strand', 'TE_class', 'TE_family', 'perc_div', 'other_column']
    df = pd.read_csv(bed_file, delim_whitespace=True, names=columns)

    df['insertion_length'] = df['length']
    df['TE_class'] = df['TE_class'].apply(lambda x: 'Satellite' if x == 'Simple_repeat' else x)
    df['TE_class'] = df['TE_class'].apply(lambda x: 'Unknown' if not isinstance(x, str) else x)
    
    return df

# Function to extract genome size from the summary file
def extract_genome_size(summary_file):
    if summary_file.endswith('.gz'):
        with gzip.open(summary_file, 'rt') as f:
            summary_lines = f.readlines()
    else:
        with open(summary_file, 'r') as f:
            summary_lines = f.readlines()

    for line in summary_lines:
        match = re.search(r'Total Length:\s+(\d+)', line)
        if match:
            genome_size = int(match.group(1))
            return genome_size
    raise ValueError("Genome size (Total Length) not found in the summary file.")

# Function to process either file type (RepeatMasker .out.gz or .bed file)
def process_file(file_type, input_file):
    if file_type == 'out':
        return parse_repeatmasker_outfile(input_file)
    elif file_type == 'bed':
        return parse_bed_file(input_file)
    else:
        raise ValueError(f"Unsupported file type: {file_type}")

# Parallelized function to calculate TE proportions per group
def calculate_group(group, genome_size):
    total_insertion = group['insertion_length'].sum()
    proportion = (total_insertion / genome_size) * 100  # Multiply by 100 to convert to percentage
    return group['binned_div'].iloc[0], group['TE_class'].iloc[0], proportion

# Function to calculate the proportion of genome occupied by TE classes
def calculate_te_proportions(df, genome_size, bin_size, selected_classes, mutation_rate=None, num_procs=1):
    if mutation_rate:
        df['binned_div'] = df['perc_div'] / mutation_rate
    else:
        df['binned_div'] = (df['perc_div'] // bin_size) * bin_size
    
    df = df[df['TE_class'].isin(selected_classes)]
    grouped = df.groupby(['binned_div', 'TE_class'])
    
    if num_procs > 1:
        with mp.Pool(processes=num_procs) as pool:
            results = pool.starmap(calculate_group, [(group, genome_size) for name, group in grouped])
    else:
        results = [calculate_group(group, genome_size) for name, group in grouped]

    te_groups = pd.DataFrame(results, columns=['binned_div', 'TE_class', 'proportion'])
    return te_groups

# Function to create stacked bar plot
# Function to create stacked bar plot and line plot
def plot_te_proportions(te_groups, output_file, bin_size, max_divergence, spacing, mutation_rate=None):
    pivot_data = te_groups.pivot(index='binned_div', columns='TE_class', values='proportion').fillna(0)
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Generate stacked bar plot
    color_list = [class_colors.get(te_class, '#000000') for te_class in pivot_data.columns]
    pivot_data.plot(kind='bar', stacked=True, width=spacing, ax=ax, color=color_list)

    # Add line plot for each TE class
    for te_class in pivot_data.columns:
        ax.plot(pivot_data.index, pivot_data[te_class], marker='o', color=class_colors.get(te_class, '#000000'), label=te_class)

    # Set labels, title, and limits
    x_label = 'Time (millions of years)' if mutation_rate else f'Genetic Divergence (%) (Bins of {bin_size})'
    ax.set_xlabel(x_label)
    ax.set_ylabel('Proportion of Genome Occupied')
    ax.set_title('Repetitive Proportion by Genetic Divergence/Class')
    ax.set_xlim([-1, max_divergence // bin_size + 1])  # Adjust x-axis limit based on max_divergence
    
    # Show legend and save the plot
    ax.legend(loc='upper right', bbox_to_anchor=(1.15, 1))
    plt.tight_layout()
    plt.savefig(output_file)

# Add this function to save the table
def save_te_proportions_table(te_groups, output_table_file):
    te_groups.to_csv(output_table_file, index=False)
    print(f"Saved TE proportions table to {output_table_file}")

# Function to create pie chart for TE classes
def plot_te_pie_chart(df, genome_size, output_file):
    te_class_proportions = df.groupby('TE_class')['insertion_length'].sum() / genome_size
    te_occupied = te_class_proportions.sum()
    non_te_proportion = 1.0 - te_occupied
    te_class_proportions = pd.concat([te_class_proportions, pd.Series([non_te_proportion], index=['Non-TE'])])

    plt.figure(figsize=(8, 8))
    color_list = [class_colors.get(te_class, '#000000') for te_class in te_class_proportions.index]
    plt.pie(te_class_proportions, labels=te_class_proportions.index, autopct='%1.1f%%', startangle=90, colors=color_list)
    plt.title('Proportion of Genome Occupied by TE Classes and Non-TE Region')
    plt.tight_layout()
    plt.savefig(output_file)

# Function to process and plot results based on file type
def process_and_plot(file_type, input_file, genome_size, output_file_prefix, bin_size, max_divergence, spacing, mutation_rate, selected_classes, num_procs, plot_type):
    df = process_file(file_type, input_file)
    
    if plot_type in ['bar', 'both']:
        te_proportions = calculate_te_proportions(df, genome_size, bin_size, selected_classes, mutation_rate, num_procs)
        plot_te_proportions(te_proportions, f"{output_file_prefix}.stackedbar.png", bin_size, max_divergence, spacing, mutation_rate)
        save_te_proportions_table(te_proportions, f"{output_file_prefix}.table.csv")
    
    if plot_type in ['pie', 'both']:
        plot_te_pie_chart(df, genome_size, f"{output_file_prefix}.pie.png")

def main():
    parser = argparse.ArgumentParser(description='Generate stacked bar and pie charts of TE proportions from RepeatMasker .out or .bed files.')
    parser.add_argument('-f', '--file_type', required=True, choices=['out', 'bed'], help='Type of file to process (out or bed)')
    parser.add_argument('-i', '--input_file', required=True, help='Input file (.out.gz or .bed)')
    parser.add_argument('-s', '--summary', required=True, help='Summary file with genome size (supports .gz)')
    parser.add_argument('-o', '--output', required=True, help='Output prefix for generated files')
    parser.add_argument('--plot_type', default='both', choices=['bar', 'pie', 'both'], help='Type of plot(s) to generate')
    parser.add_argument('-bin', '--bin_size', type=int, default=1, help='Bin size for genetic divergence')
    parser.add_argument('-m', '--max_divergence', type=int, default=50, help='Maximum divergence for x-axis')
    parser.add_argument('--spacing', type=float, default=0.8, help='Bar width for stacked bar chart')
    parser.add_argument('-mr', '--mutation_rate', type=float, help='Mutation rate to estimate time since divergence')
    parser.add_argument('-c', '--selected_classes', nargs='+', default=['DNA', 'LINE', 'SINE', 'LTR', 'RC', 'Satellite', 'Unknown'], help='List of TE classes to include in the plot')
    #Example alternative to default: --selected_classes DNA,DIRS,LINE,LTR,RC,SINE,NonLTR,Unknown,Other
    parser.add_argument('--num_procs', type=int, default=1, help='Number of processors for parallel processing')

    args = parser.parse_args()

    genome_size = extract_genome_size(args.summary)
    process_and_plot(
        file_type=args.file_type,
        input_file=args.input_file,
        genome_size=genome_size,
        output_file_prefix=args.output,
        bin_size=args.bin_size,
        max_divergence=args.max_divergence,
        spacing=args.spacing,
        mutation_rate=args.mutation_rate,
        selected_classes=args.selected_classes,
        num_procs=args.num_procs,
        plot_type=args.plot_type
    )

if __name__ == '__main__':
    main()

