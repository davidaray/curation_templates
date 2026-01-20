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
    te_class_proportions = pd.concat([te_class_proportions, pd.Series([non_te_proportion], index=['Non-TE'])])

    plt.figure(figsize=(8, 8))
    color_list = [class_colors.get(te_class, '#000000') for te_class in te_class_proportions.index]
    plt.pie(te_class_proportions, labels=te_class_proportions.index, autopct='%1.1f%%', startangle=90, colors=color_list)
    plt.title('Proportion of Genome Occupied by TE Classes and Non-TE Region')
    plt.tight_layout()
    plt.savefig(output_file)
    
def process_and_plot(df, genome_size, output_file_prefix, bin_size, max_divergence, spacing, mutation_rate, selected_classes, num_procs, plot_type):
    if plot_type in ['bar', 'both']:
        te_proportions = calculate_te_proportions(df, genome_size, bin_size, selected_classes, mutation_rate, num_procs)
        plot_te_proportions(te_proportions, f"{output_file_prefix}.stackedbar.png", bin_size, max_divergence, spacing, mutation_rate)
        save_te_proportions_table(te_proportions, f"{output_file_prefix}.table.csv")  # Save table with the TE proportions
    
    if plot_type in ['pie', 'both']:
        plot_te_pie_chart(df, genome_size, f"{output_file_prefix}.pie.png")

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


'''
Use argparse in python to take in a repeatmasker .out file (-r) and generate a stacked bar plot (-o) of the data. 

Example data from .out file:
   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat
score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID

  371  34.0  5.9  5.9  Chr19        2682    3037  (156962) +  pMex_TIR_123-  DNA/Unknown            211  566  (988)      1
  670  33.2  4.6  4.6  Chr19        9190    9858  (150141) +  L1-1_PTr       LINE/L1               3584 4252 (1843)      2
  732  33.6  4.3  2.1  Chr19       10766   11227  (148772) +  pMex_TIR_123-  DNA/Unknown            199  670  (884)      3
  554  33.8  4.8  4.8  Chr19       15193   15861  (144138) +  L1-1_PTr       LINE/L1               3584 4252 (1843)      4
  257  21.3  0.0  0.0  Chr19       16415   16461  (143538) C  L1-1_PTr_pMex-5 LINE/L1             (4887) 1278   1232      5
  850  33.3  2.4  2.6  Chr19       16783   17241  (142758) +  pMex_TIR_123-  DNA/Unknown            213  670  (884)      6
  232  25.0  0.0  2.9  Chr19       17823   17892  (142107) +  pMex_TIR_123-  DNA/Unknown            866  933  (621)      7
   12   8.9  4.2  0.0  Chr19       19666   19689  (140310) +  (CTT)n         Simple_repeat            1   25    (0)      8
  256  12.5  4.5  4.5  Chr19       21433   21499  (138500) +  pMex_TIR_123-  DNA/Unknown           1335 1401  (153)      9
  512  28.6  3.3  3.3  Chr19       21985   22259  (137740) C  L1-1_PTr_pMex-5 LINE/L1             (3677) 2488   2214     10
  274  36.9  1.5  0.0  Chr19       22531   22660  (137339) +  L1-1_PTr       LINE/L1               4140 4271 (1824)     11
  337  35.4  4.1  4.1  Chr19       23562   23906  (136093) +  pMex_TIR_123-  DNA/Unknown            211  555  (999)     12
  548  32.8  4.4  4.7  Chr19       28060   28716  (131283) +  L1-1_PTr       LINE/L1               3600 4254 (1841)     13
 2094  28.6  2.9  2.1  Chr19       29423   30208  (129791) +  pMex_TIR_123-  DNA/Unknown              5  796  (758)     14

The x-axis should display genetic divergence, column 2 of the .out file. The y-axis should display proportionof the genome occupied by each class of transposable element (TE) in the table. The number of bases occupied by each TE insertion can be found by calculating the difference between columns 6 and 7. The class of each TE is the text before '/' in column 11.

Genome size in base paris (-g) should be input using argparse. To calculate proportion of the genome, the sum of all TE insertion base pairs divided by the genome size should be determined.

Use your best judgement to determine bin sizes on the x-axis and height of the y-axis. 

Allow x-axis bin sizes to be set via argparse.

Allow the .out file to be zipped as a .gz file. Keep the same input option value. Also grab input from the summary.gz file provided (-s). This first several lines of the summary file will have this structure: 
Repeat Classes
==============
Total Sequences: 1
Total Length: 159999 bp
Class                  Count        bpMasked    %masked
=====                  =====        ========     =======
The genome size is equal to the Total Length. In this case, 15999.

Starting with the script below, add the following. 
1. The option to choose either a pie chart or stacked bar plot or both via argparse.
2. If the pie chart option is chosen, generate a pie chart using the TE class data and the total genome proportions of each class relative to the genome as a whole. Make region of the pie chart not occupied by transposable elements black.

Please incorporate the option for a neutral mutation rate (-n).

If the neutral mutation rate is provided, it should be in the format of xe-y, meaning x times 10 to the -y power. 

If the mutation rate is provided, calculate values for time on the x-axis using the equation, genetic divergence = mutation rate x time. 

If the mutation rate is not provided, use genetic divergence as before. 

Add an option to define the spacing of the bins (-sp). Right now, they're too far apart. Also make the default maximum genetic divergence and the corresponding time on the x-axis, 40 with the option to change it (-maxd).

Keep everything as is but print the percent, non-TE fraction in the black portion similar to the way the TE portions are printed on the pie chart.

Parallelize for memory usage and speed improvement, please.

1. Allow the user to indicate what classes should be included in the plots using an option --classes. The input format for classes should be text as shown in the following example (DNA, DIRS, LINE, LTR, RC, SINE, NonLTR, Satellite, Unknown, Other'). These would be the default categories.  Satellite and Simple Repeat should be combined into the 'Satellite' category. 'Other' would include RNA, Segmenta, rRNA, scRNA, snRNA, and tRNA. 'NonLTR' would include NonLTR, Retrogene, and Retrotransposon.

2. Change the stacked bar graph title to "Repetitive Proportion by Genetic Divergence/Class".

Now, save the pie chart and stacked bar plot as separate files. Currently, the pie chart is overwriting the stacked bar. Incorporate 'pie' and 'stackedbar' into the output file name provided by the user.

Remove the lines that Center the percentage fraction of Non-TE in the black slice.

Need to incorporate:
The bins when using the mutation rate to calculate time are causing problems when -n is included. Convert the bin size as appropriate to match the bin sizes displayed when -n is not included. 

'''
