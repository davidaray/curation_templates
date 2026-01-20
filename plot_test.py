import argparse
import pandas as pd
import matplotlib.pyplot as plt
import gzip
import re
import multiprocessing as mp

# Use non-interactive Agg backend
import matplotlib
matplotlib.use('Agg')

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
    
    def categorize_te_class(te_class):
        if te_class in ['Satellite', 'Simple_repeat']:
            return 'Satellite'
        elif 'DNA' in te_class:
            return 'DNA'
        elif 'LINE' in te_class:
            return 'LINE'
        elif 'LTR' in te_class:
            return 'LTR'
        elif 'SINE' in te_class:
            return 'SINE'
        elif 'NonLTR' in te_class or 'Retrogene' in te_class or 'Retroposon' in te_class or 'Retrotransposon' in te_class:
            return 'NonLTR'
        elif 'RC' in te_class:
            return 'RC'
        elif 'RNA' in te_class or 'rRNA' in te_class or 'scRNA' in te_class or 'snRNA' in te_class or 'tRNA' in te_class:
            return 'RNA'
        elif te_class == 'Unknown':
            return 'Unknown'
        else:
            return 'Other'

    df['TE_class'] = df['TE_class'].apply(categorize_te_class)
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
def plot_te_proportions(te_groups, output_file, bin_size, max_divergence, spacing, mutation_rate=None):
    pivot_data = te_groups.pivot(index='binned_div', columns='TE_class', values='proportion').fillna(0)
    fig, ax = plt.subplots(figsize=(10, 6))
    color_list = [class_colors.get(te_class, '#000000') for te_class in pivot_data.columns]
    pivot_data.plot(kind='bar', stacked=True, width=spacing, ax=ax, color=color_list)

    x_label = 'Time (millions of years)' if mutation_rate else f'Genetic Divergence (%) (Bins of {bin_size})'
    ax.set_xlabel(x_label)
    ax.set_ylabel('Proportion of Genome Occupied')
    ax.set_title('Repetitive Proportion by Genetic Divergence/Class')
    ax.set_xlim([-1, max_divergence // bin_size + 1])  # Adjust x-axis limit based on max_divergence
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
    te_class_proportions = te_class_proportions.append(pd.Series([non_te_proportion], index=['Non-TE']))

    plt.figure(figsize=(8, 8))
    color_list = [class_colors.get(te_class, '#000000') for te_class in te_class_proportions.index]
    plt.pie(te_class_proportions, labels=te_class_proportions.index, autopct='%1.1f%%', startangle=90, colors=color_list)
    plt.title('Proportion of Genome Occupied by TE Classes and Non-TE Region')
    plt.tight_layout()
    plt.savefig(output_file)

def process_and_plot(df, genome_size, output_file_prefix, bin_size, max_divergence, spacing, mutation_rate, selected_classes, num_procs, plot_type):
    if plot_type in ['bar', 'both']:
        te_proportions = calculate_te_proportions(df, genome_size, bin_size, selected_classes, mutation_rate, num_procs)
        plot_te_proportions(te_proportions, f"{output_file_prefix}_stackedbar.png", bin_size, max_divergence, spacing, mutation_rate)
        save_te_proportions_table(te_proportions, f"{output_file_prefix}_table.csv")  # Save table with the TE proportions
    
    if plot_type in ['pie', 'both']:
        plot_te_pie_chart(df, genome_size, f"{output_file_prefix}_pie.png")

def main():
    parser = argparse.ArgumentParser(description='Generate a stacked bar and pie chart of TE proportions.')
    parser.add_argument('-r', '--repeatmasker', required=True, help='RepeatMasker .out file (supports .gz)')
    parser.add_argument('-b', '--bed', required=True, help='RepeatMasker .bed file')
    parser.add_argument('-s', '--summary', required=True, help='Summary file with genome size (supports .gz)')
    parser.add_argument('-o', '--output', required=True, help='Output prefix for generated files')
    parser.add_argument('--plot_type', default='both', choices=['bar', 'pie', 'both'], help='Type of plot(s) to generate')
    parser.add_argument('-bin', '--bin_size', type=int, default=1, help='Bin size for genetic divergence')
    parser.add_argument('--max_divergence', type=int, default=50, help='Maximum divergence for the x-axis')
    parser.add_argument('--spacing', type=float, default=0.8, help='Spacing for bars')
    parser.add_argument('--mutation_rate', type=float, help='Mutation rate for calculating time')
    parser.add_argument('--classes', default='all', help='Comma-separated list of TE classes to include (default: all)')
    #Example alternative to all: --classes DNA,DIRS,LINE,LTR,RC,SINE,NonLTR,Unknown,Other
    parser.add_argument('-proc', '--num_procs', type=int, default=1, help='Number of processors to use')
    args = parser.parse_args()

    selected_classes = args.classes.split(',') if args.classes != 'all' else list(class_colors.keys())

    genome_size = extract_genome_size(args.summary)
    
    # Process .out file
    df_out = parse_repeatmasker_outfile(args.repeatmasker)
    process_and_plot(df_out, genome_size, f"{args.output}_out", args.bin_size, args.max_divergence, args.spacing, args.mutation_rate, selected_classes, args.num_procs, args.plot_type)
    
    # Process .bed file
    df_bed = parse_bed_file(args.bed)
    process_and_plot(df_bed, genome_size, f"{args.output}_bed", args.bin_size, args.max_divergence, args.spacing, args.mutation_rate, selected_classes, args.num_procs, args.plot_type)

if __name__ == '__main__':
    main()
