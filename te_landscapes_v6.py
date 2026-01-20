Change this script to add an argument on where to find the .bed files. 

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from collections import defaultdict

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Analyze transposable element landscapes by divergence'
    )
    parser.add_argument(
        '-g', '--genomesize',
        required=True,
        help='Genome size file with format: ID, genome_size_bp, mutation_rate, Species_ID'
    )
    parser.add_argument(
        '-d', '--divergence',
        type=float,
        default=50.0,
        help='Maximum divergence threshold (default: 50)'
    )
    parser.add_argument(
        '-m', '--minimum100bp',
        action='store_true',
        help='Include hits under 100bp (default: exclude them)'
    )
    
    return parser.parse_args()

def read_genome_sizes(filename):
    """Read genome size file and return dictionary of genome info"""
    genomes = {}
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 4:
                        genome_id = parts[0]
                        genome_size = int(parts[1])
                        mutation_rate = float(parts[2])
                        species_id = parts[3]
                        genomes[genome_id] = {
                            'size': genome_size,
                            'mutation_rate': mutation_rate,
                            'species': species_id
                        }
    except FileNotFoundError:
        print(f"Error: Genome size file '{filename}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading genome size file: {e}")
        sys.exit(1)
    
    return genomes

def classify_te_type(te_name, class_col, family_col):
    """
    Classify transposable element based on name and class/family information
    Returns one of: DNA, RC, LINE, SINE, LTR, Unknown
    """
    te_name = str(te_name).upper()
    class_col = str(class_col).upper()
    family_col = str(family_col).upper()
    
    # Simple repeats and satellites - exclude from landscape
    if 'SIMPLE_REPEAT' in class_col or ')N' in te_name or 'SATELLITE' in class_col:
        return None
    
    # LTR retrotransposons
    if 'LTR' in class_col or 'LTR' in te_name:
        return 'LTR'
    
    # LINE elements
    if 'LINE' in class_col or 'NONLTR' in class_col or 'NON-LTR' in class_col:
        return 'LINE'
    
    # SINE elements
    if 'SINE' in class_col:
        return 'SINE'
    
    # DNA transposons
    if 'DNA' in class_col or 'TIR' in te_name:
        return 'DNA'
    
    # Rolling circle (RC) transposons
    if 'RC' in class_col or 'HELITRON' in class_col:
        return 'RC'
    
    # Default to Unknown
    return 'Unknown'

def process_bed_file(filename, genome_id, genome_size, max_divergence, include_small):
    """Process a single BED file and return TE landscape data"""
    
    try:
        # Try to read the BED file
        df = pd.read_csv(filename, sep='\t', header=None)
        
        # Expected columns: Scaffold, Start, End, TE, Hit_size, Orientation, Class, Family, Divergence, RM_ID
        if len(df.columns) < 10:
            print(f"Warning: {filename} has fewer than 10 columns. Skipping.")
            return None
            
        df.columns = ['Scaffold', 'Start', 'End', 'TE', 'Hit_size', 'Orientation', 'Class', 'Family', 'Divergence', 'RM_ID']
        
        # Convert divergence to numeric, handling non-numeric values
        df['Divergence'] = pd.to_numeric(df['Divergence'], errors='coerce')
        
        # Filter by divergence (exclude negative values like -1.0 for simple repeats)
        df = df[(df['Divergence'] >= 0) & (df['Divergence'] <= max_divergence)]
        
        # Filter by size if needed
        if not include_small:
            df = df[df['Hit_size'] >= 100]
        
        # Create divergence bins (0-1%, 1-2%, etc.)
        max_bin = int(max_divergence) + 1
        bins = list(range(0, max_bin + 1))
        df['Divergence_Bin'] = pd.cut(df['Divergence'], bins=bins, right=False, include_lowest=True)
        
        # Initialize landscape data structure
        landscape_data = defaultdict(lambda: defaultdict(int))
        
        # Process each row
        for _, row in df.iterrows():
            te_type = classify_te_type(row['TE'], row['Class'], row['Family'])
            
            if te_type:  # Only process if it's a valid TE type
                bin_label = row['Divergence_Bin']
                if pd.notna(bin_label):
                    bin_start = int(bin_label.left)
                    landscape_data[te_type][bin_start] += row['Hit_size']
        
        # Convert to proportions
        landscape_proportions = {}
        for te_type in landscape_data:
            landscape_proportions[te_type] = {}
            for bin_start in landscape_data[te_type]:
                proportion = (landscape_data[te_type][bin_start] / genome_size) * 100
                landscape_proportions[te_type][bin_start] = proportion
        
        return landscape_proportions
        
    except FileNotFoundError:
        print(f"Warning: BED file '{filename}' not found for genome {genome_id}")
        return None
    except Exception as e:
        print(f"Error processing {filename}: {e}")
        return None

def create_landscape_plot(landscape_data, genome_id, ax, max_divergence):
    """Create a landscape line plot for a single genome"""
    
    # Define colors for TE classes
    te_colors = {
        'DNA': '#FF0000',      # red
        'RC': '#800080',       # purple
        'LINE': '#0000FF',     # blue
        'SINE': '#ADD8E6',     # light blue
        'LTR': '#008000',      # green
        'Unknown': '#2F2F2F'   # dark grey
    }
    
    # Define the order for plotting
    te_order = ['DNA', 'RC', 'LINE', 'SINE', 'LTR', 'Unknown']
    
    # Prepare data for plotting
    max_bin = int(max_divergence)
    x_bins = list(range(0, max_bin))
    
    # Initialize data matrix
    plot_data = {}
    for te_type in te_order:
        plot_data[te_type] = [landscape_data.get(te_type, {}).get(bin_val, 0) for bin_val in x_bins]
    
    # Create line plots for each TE type
    for te_type in te_order:
        if te_type in landscape_data and any(plot_data[te_type]):
            ax.plot(x_bins, plot_data[te_type], 
                   color=te_colors[te_type], linewidth=2, label=te_type)
    
    # Formatting
    ax.set_xlabel('Divergence (%)', fontsize=10)
    ax.set_ylabel('Genome Proportion (%)', fontsize=10)
    ax.set_title(f'{genome_id}', fontsize=12, fontweight='bold')
    ax.set_xlim(0, max_bin)
    ax.grid(True, alpha=0.3)
    
    # Set x-axis ticks every 5%
    x_ticks = list(range(0, max_bin + 1, 5))
    ax.set_xticks(x_ticks)
    
    # Add legend to each individual plot
    ax.legend(loc='upper right', fontsize=8)
    
    return plot_data

def create_individual_landscape(landscape_data, genome_id, output_prefix, max_divergence):
    """Create and save an individual landscape plot"""
    
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_data = create_landscape_plot(landscape_data, genome_id, ax, max_divergence)
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_{genome_id}_landscape.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_{genome_id}_landscape.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return plot_data

def create_combined_landscapes(all_results, output_prefix, max_divergence):
    """Create a single multipanel figure with all landscape plots"""
    
    n_samples = len(all_results)
    
    # Calculate subplot layout
    cols = min(3, n_samples)
    rows = int(np.ceil(n_samples / cols))
    
    # Create figure
    fig_width = cols * 6 + 3  # Extra space for legend
    fig_height = rows * 5
    fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(fig_width, fig_height))
    
    # Handle single subplot case
    if n_samples == 1:
        axes = [axes]
    elif rows == 1:
        axes = axes.flatten()
    else:
        axes = axes.flatten()
    
    # Define colors for legend
    te_colors = {
        'DNA': '#FF0000',      # red
        'RC': '#800080',       # purple
        'LINE': '#0000FF',     # blue
        'SINE': '#ADD8E6',     # light blue
        'LTR': '#008000',      # green
        'Unknown': '#2F2F2F'   # dark grey
    }
    
    # Create landscape plots for each sample (without individual legends)
    all_plot_data = {}
    
    for i, (sample_id, landscape_data) in enumerate(all_results.items()):
        ax = axes[i]
        
        # Define the order for plotting
        te_order = ['DNA', 'RC', 'LINE', 'SINE', 'LTR', 'Unknown']
        
        # Prepare data for plotting
        max_bin = int(max_divergence)
        x_bins = list(range(0, max_bin))
        
        # Initialize data matrix
        plot_data = {}
        for te_type in te_order:
            plot_data[te_type] = [landscape_data.get(te_type, {}).get(bin_val, 0) for bin_val in x_bins]
        
        # Create line plots for each TE type
        for te_type in te_order:
            if te_type in landscape_data and any(plot_data[te_type]):
                ax.plot(x_bins, plot_data[te_type], 
                       color=te_colors[te_type], linewidth=2, label=te_type)
        
        # Formatting
        ax.set_xlabel('Divergence (%)', fontsize=10)
        ax.set_ylabel('Genome Proportion (%)', fontsize=10)
        ax.set_title(f'{sample_id}', fontsize=12, fontweight='bold')
        ax.set_xlim(0, max_bin)
        ax.grid(True, alpha=0.3)
        
        # Set x-axis ticks every 5%
        x_ticks = list(range(0, max_bin + 1, 5))
        ax.set_xticks(x_ticks)
        
        all_plot_data[sample_id] = plot_data
    
    # Hide unused subplots
    for j in range(n_samples, len(axes)):
        axes[j].set_visible(False)
    
    # Create legend with all TE types that appear in the data
    legend_elements = []
    legend_labels = []
    
    te_order = ['DNA', 'RC', 'LINE', 'SINE', 'LTR', 'Unknown']
    for te_type in te_order:
        # Check if this TE type appears in any sample
        appears = False
        for sample_data in all_results.values():
            if te_type in sample_data and any(sample_data[te_type].values()):
                appears = True
                break
        
        if appears:
            legend_elements.append(plt.Line2D([0], [0], color=te_colors[te_type], linewidth=2))
            legend_labels.append(te_type)
    
    fig.legend(legend_elements, legend_labels, loc='center right', 
              bbox_to_anchor=(0.98, 0.5), fontsize=12)
    
    plt.tight_layout()
    plt.subplots_adjust(right=0.85)  # Make room for legend
    
    # Save the figure
    plt.savefig(f'{output_prefix}_te_landscapes.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_te_landscapes.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    return all_plot_data

def save_landscape_tables(all_plot_data, output_prefix, max_divergence):
    """Save landscape data tables for all samples"""
    
    max_bin = int(max_divergence)
    x_bins = list(range(0, max_bin))
    
    # Create summary table
    for sample_id, plot_data in all_plot_data.items():
        # Create individual sample table
        table_data = {'Divergence_Bin': [f'{i}-{i+1}%' for i in x_bins]}
        
        te_order = ['DNA', 'RC', 'LINE', 'SINE', 'LTR', 'Unknown']
        for te_type in te_order:
            if te_type in plot_data:
                table_data[te_type] = [f'{val:.4f}' for val in plot_data[te_type]]
            else:
                table_data[te_type] = ['0.0000'] * len(x_bins)
        
        # Save individual table
        df = pd.DataFrame(table_data)
        df.to_csv(f'{output_prefix}_{sample_id}_landscape.csv', index=False)
    
    # Create combined summary table
    all_samples_data = []
    for sample_id, plot_data in all_plot_data.items():
        for i, bin_start in enumerate(x_bins):
            row = {
                'Sample_ID': sample_id,
                'Divergence_Bin': f'{bin_start}-{bin_start+1}%',
                'Bin_Start': bin_start
            }
            
            te_order = ['DNA', 'RC', 'LINE', 'SINE', 'LTR', 'Unknown']
            for te_type in te_order:
                if te_type in plot_data:
                    row[te_type] = plot_data[te_type][i]
                else:
                    row[te_type] = 0.0
            
            all_samples_data.append(row)
    
    combined_df = pd.DataFrame(all_samples_data)
    combined_df.to_csv(f'{output_prefix}_all_samples_landscape.csv', index=False)
    
    print(f"Landscape tables saved:")
    for sample_id in all_plot_data.keys():
        print(f"- {output_prefix}_{sample_id}_landscape.csv")
    print(f"- {output_prefix}_all_samples_landscape.csv")

def main():
    args = parse_arguments()
    
    # Read genome size information
    genomes = read_genome_sizes(args.genomesize)
    
    if not genomes:
        print("No genomes found in genome size file.")
        sys.exit(1)
    
    print(f"Found {len(genomes)} genomes to process")
    print(f"Maximum divergence: {args.divergence}")
    print(f"Include hits <100bp: {args.minimum100bp}")
    
    # Process each genome
    all_results = {}
    all_plot_data = {}
    successful_genomes = []
    
    output_prefix = "te_landscape"
    
    for genome_id, genome_info in genomes.items():
        bed_filename = f"{genome_id}_rm.bed"
        print(f"Processing {genome_id} ({bed_filename})...")
        
        result = process_bed_file(
            bed_filename, 
            genome_id, 
            genome_info['size'], 
            args.divergence, 
            args.minimum100bp
        )
        
        if result:
            all_results[genome_id] = result
            successful_genomes.append(genome_id)
            
            # Create individual landscape plot
            print(f"  Creating individual landscape plot for {genome_id}...")
            plot_data = create_individual_landscape(result, genome_id, output_prefix, args.divergence)
            all_plot_data[genome_id] = plot_data
    
    if not all_results:
        print("No data processed successfully.")
        sys.exit(1)
    
    # Create multipanel landscape plot
    print(f"\nCreating combined landscape plot with {len(all_results)} samples...")
    create_combined_landscapes(all_results, output_prefix, args.divergence)
    
    # Save landscape tables
    print("Saving landscape data tables...")
    save_landscape_tables(all_plot_data, output_prefix, args.divergence)
    
    print(f"\nOutput files created:")
    print("- te_landscape_te_landscapes.pdf")
    print("- te_landscape_te_landscapes.png")
    for genome_id in successful_genomes:
        print(f"- te_landscape_{genome_id}_landscape.pdf")
        print(f"- te_landscape_{genome_id}_landscape.png")
    print("- Individual and combined CSV tables")
    
    print(f"\nAnalysis complete! Processed {len(all_results)} samples.")

if __name__ == "__main__":
    main()

