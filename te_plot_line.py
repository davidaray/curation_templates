import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import gzip
import re
import multiprocessing as mp

class_colors = {
    'DNA': '#1f77b4', 'DIRS': '#ff7f0e', 'LINE': '#2ca02c', 'LTR': '#d62728', 
    'RC': '#9467bd', 'SINE': '#8c564b', 'NonLTR': '#e377c2', 'Satellite': '#7f7f7f', 
    'Unknown': '#bcbd22', 'Other': '#17becf', 'Non-TE': '#000000'
}

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

def parse_bed_file(bed_file):
    columns = ['scaffold', 'start', 'end', 'repeat_name', 'length', 'strand', 'TE_class', 'TE_family', 'perc_div', 'other_column']
    df = pd.read_csv(bed_file, delim_whitespace=True, names=columns)
    df['insertion_length'] = df['length']
    df['TE_class'] = df['TE_class'].apply(lambda x: 'Satellite' if x == 'Simple_repeat' else x)
    df['TE_class'] = df['TE_class'].apply(lambda x: 'Unknown' if not isinstance(x, str) else x)
    return df

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

def calculate_group(group, genome_size):
    total_insertion = group['insertion_length'].sum()
    proportion = (total_insertion / genome_size) * 100
    return group['binned_div'].iloc[0], group['TE_class'].iloc[0], proportion

def calculate_te_proportions(df, genome_size, bin_size, selected_classes, mutation_rate=None, num_procs=1):
    if mutation_rate:
        df['binned_div'] = df['perc_div'] / mutation_rate
    else:
        df['binned_div'] = (df['perc_div'] // bin_size) * bin_size

    df['TE_class'] = df['TE_class'].apply(lambda x: x if x in selected_classes else 'Other')
    grouped = df.groupby(['binned_div', 'TE_class'])
    
    if num_procs > 1:
        with mp.Pool(processes=num_procs) as pool:
            results = pool.starmap(calculate_group, [(group, genome_size) for name, group in grouped])
    else:
        results = [calculate_group(group, genome_size) for name, group in grouped]

    te_groups = pd.DataFrame(results, columns=['binned_div', 'TE_class', 'proportion'])
    return te_groups

def plot_te_proportions(te_groups, output_file, bin_size, max_divergence, mutation_rate=None):
    pivot_data = te_groups.pivot(index='binned_div', columns='TE_class', values='proportion').fillna(0)
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for te_class in pivot_data.columns:
        ax.plot(pivot_data.index, pivot_data[te_class], label=te_class, color=class_colors.get(te_class, '#000000'))

    x_label = 'Time (millions of years)' if mutation_rate else f'Genetic Divergence (%) (Bins of {bin_size})'
    ax.set_xlabel(x_label)
    ax.set_ylabel('Proportion of Genome Occupied')
    ax.set_title('Repetitive Proportion by Genetic Divergence/Class')
    ax.set_xlim([-1, max_divergence // bin_size + 1])
    ax.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(output_file)

def save_te_proportions_table(te_groups, output_table_file):
    te_groups.to_csv(output_table_file, index=False)
    print(f"Saved TE proportions table to {output_table_file}")

def process_and_plot(df, genome_size, output_file_prefix, bin_size, max_divergence, mutation_rate, selected_classes, num_procs):
    te_proportions = calculate_te_proportions(df, genome_size, bin_size, selected_classes, mutation_rate, num_procs)
    plot_te_proportions(te_proportions, f"{output_file_prefix}.lineplot.png", bin_size, max_divergence, mutation_rate)
    save_te_proportions_table(te_proportions, f"{output_file_prefix}.table.csv")

def main():
    parser = argparse.ArgumentParser(description='Generate a line plot of TE proportions.')
    parser.add_argument('-r', '--repeatmasker', required=True, help='RepeatMasker .out file (supports .gz)')
    parser.add_argument('-b', '--bed', required=True, help='RepeatMasker .bed file')
    parser.add_argument('-s', '--summary', required=True, help='Summary file with genome size (supports .gz)')
    parser.add_argument('-o', '--output', required=True, help='Output prefix for generated files')
    parser.add_argument('-bin', '--bin_size', type=int, default=1, help='Bin size for genetic divergence')
    parser.add_argument('--max_divergence', type=int, default=50, help='Maximum divergence for the x-axis')
    parser.add_argument('--mutation_rate', type=float, help='Mutation rate for calculating time')
    parser.add_argument('--classes', default='all', help='Comma-separated list of TE classes to include (default: all)')
    parser.add_argument('-proc', '--num_procs', type=int, default=1, help='Number of processors to use')
    args = parser.parse_args()

    selected_classes = args.classes.split(',') if args.classes != 'all' else list(class_colors.keys())
    genome_size = extract_genome_size(args.summary)
    
    df_out = parse_repeatmasker_outfile(args.repeatmasker)
    process_and_plot(df_out, genome_size, f"{args.output}_out", args.bin_size, args.max_divergence, args.mutation_rate, selected_classes, args.num_procs)
    
    df_bed = parse_bed_file(args.bed)
    process_and_plot(df_bed, genome_size, f"{args.output}_bed", args.bin_size, args.max_divergence, args.mutation_rate, selected_classes, args.num_procs)

if __name__ == '__main__':
    main()

""" 
Changes Made when compared to te_plot_stacked_pie.py
Replaced Stacked Bar Plot with Line Plot: The plot_te_proportions function now generates a line plot, where each line represents a TE class, colored based on class_colors. Removed 
Pie Chart: The plot_te_pie_chart function and any references to it were removed.
Adjusted Output File Naming: The line plot output files are saved as output_prefix.lineplot.png.
"""
