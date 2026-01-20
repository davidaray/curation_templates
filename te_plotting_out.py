import argparse
import gzip
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from collections import defaultdict
import numpy as np
import concurrent.futures

# Define default TE classes and colors
TE_CLASSES = ["LINE", "SINE", "LTR", "DIRS", "DNA", "RC", "Unknown", "Satellite", "Simple_repeat", "Other", "NonLTR"]
COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8"]

def parse_args():
    parser = argparse.ArgumentParser(description="Generate TE plots from RepeatMasker .out file.")
    parser.add_argument("-io", "--input_out_file", required=True, help="Input RepeatMasker .out.gz file")
    parser.add_argument("-s", "--summary", required=True, help="Input summary.gz file containing genome size")
    parser.add_argument("-op", "--output_prefix", required=True, help="Prefix for output files")
    parser.add_argument("-bin", "--bin_size", type=int, default=1, help="Bin size for genetic divergence (default: 1)")
    parser.add_argument("-sp", "--spacing", type=float, default=0.8, help="Column widths for stacked bar plot (default: 0.8)")
    parser.add_argument("-c", "--selected_classes", default=",".join(TE_CLASSES), help="Comma-separated TE classes to include")
    parser.add_argument("-proc", "--num_proc", type=int, default=1, help="Number of processors for parallel processing")
    return parser.parse_args()

def get_genome_size(summary_file):
    with gzip.open(summary_file, 'rt') as f:
        for i, line in enumerate(f):
            if i == 3:  # Genome size is on the fourth line
                return int(line.split(":")[1].strip().split()[0])

def read_out_file(out_file, selected_classes):
    data = []
    with gzip.open(out_file, 'rt') as f:
        for line in f.readlines()[3:]:  # Skip header lines
            fields = line.strip().split()
            div, start, end, te_class = float(fields[1]), int(fields[5]), int(fields[6]), fields[10].split('/')[0]
            if te_class in selected_classes:
                data.append((te_class, div, end - start))
    return pd.DataFrame(data, columns=['TE_Class', 'Divergence', 'Length'])

def calculate_proportions(df, genome_size, bin_size, max_divergence=50):
    bins = np.arange(0, max_divergence + bin_size, bin_size)
    proportions = defaultdict(list)
    for te_class, group in df.groupby('TE_Class'):
        hist, _ = np.histogram(group['Divergence'], bins=bins, weights=group['Length'])
        proportions[te_class] = hist / genome_size
    return pd.DataFrame(proportions, index=bins[:-1])

def plot_stacked_bar(df, classes, colors, output_file, spacing):
    ax = df.plot(kind='bar', stacked=True, color=colors, width=spacing)
    ax.set_xticks(range(0, len(df.index), 10))  # Adjust interval as needed
    ax.set_xticklabels(range(0, 51, 10), fontsize=8)  # Use 0 to 50 with step 10
    plt.xlabel("Divergence")
    plt.ylabel("Genome Proportion")
    plt.title("TE Genome Proportions by Divergence (Stacked Bar)")
    plt.yticks(fontsize=8)
    plt.savefig(output_file, format='png', dpi=300)
    plt.close()

def plot_line(df, classes, colors, output_file):
    for te_class, color in zip(classes, colors):
        plt.plot(df.index, df[te_class], label=te_class, color=color)
    plt.xlabel("Divergence")
    plt.ylabel("Genome Proportion")
    plt.title("TE Genome Proportions by Divergence (Line Plot)")
    plt.legend()
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.savefig(output_file, format='png', dpi=300)
    plt.close()

def plot_pie(df, classes, colors, output_file):
    proportions = df.sum(axis=0)
    total_covered = proportions.sum()
    unmasked_proportion = 1.0 - total_covered  # Assuming proportions are in decimal form
    proportions = proportions.append(pd.Series([unmasked_proportion], index=["Unmasked"]))
    colors = colors + ["#000000"]  # Black color for "Unmasked"
    classes = classes + ["Unmasked"]
    wedges, texts, autotexts = plt.pie(
        proportions, labels=classes, colors=colors, autopct='%1.1f%%', textprops={'fontsize': 8}
    )
    for i, label in enumerate(classes):
        if label == "Unmasked":
            texts[i].set_color("white")
            autotexts[i].set_color("white")
    plt.title("TE Genome Proportions by Class (Pie Chart)", fontsize=10)
    plt.savefig(output_file, format='png', dpi=300)
    plt.close()
    
def main():
    args = parse_args()
    genome_size = get_genome_size(args.summary)
    selected_classes = args.selected_classes.split(',')
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.num_proc) as executor:
        future_out = executor.submit(read_out_file, args.input_out_file, selected_classes)
        out_df = future_out.result()
    
    out_proportions = calculate_proportions(out_df, genome_size, args.bin_size)
    out_proportions.to_csv(f"{args.output_prefix}_out_data.csv")
    plot_stacked_bar(out_proportions, selected_classes, COLORS, f"{args.output_prefix}_out_stacked_bar.png", args.spacing)
    plot_line(out_proportions, selected_classes, COLORS, f"{args.output_prefix}_out_line.png")
    plot_pie(out_proportions, selected_classes, COLORS, f"{args.output_prefix}_out_pie.png")

if __name__ == "__main__":
    main()

