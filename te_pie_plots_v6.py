#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def get_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate pie charts and tables for repeat element analysis from RepeatMasker BED files."
    )
    parser.add_argument(
        "-g", "--genomesize", 
        required=True, 
        type=str, 
        help="Path to genome size file. Tab-delimited: ID, genome_size(bp), mutation_rate, Species_ID"
    )
    parser.add_argument(
        "-d", "--divergence", 
        type=int, 
        default=50, 
        help="Maximum divergence value. Default = 50"
    )
    parser.add_argument(
        "-m", "--minimum100bp", 
        action='store_true', 
        help="If entered, count hits under 100 bp. Default is to omit them."
    )
    
    return parser.parse_args()

def resolve_label_overlaps(external_labels, min_distance=0.3, max_radius=1.6):
    """
    Resolve overlapping external labels by adjusting their positions.
    
    Args:
        external_labels: List of (angle, x1, y1, text) tuples
        min_distance: Minimum distance between label centers
        max_radius: Maximum radius for label placement
    
    Returns:
        List of (angle, x1, y1, x2, y2, text) tuples with resolved positions
    """
    if len(external_labels) <= 1:
        # No overlap possible with 0 or 1 labels
        result = []
        for angle, x1, y1, text in external_labels:
            x2 = 1.4 * np.cos(np.radians(angle))
            y2 = 1.4 * np.sin(np.radians(angle))
            result.append((angle, x1, y1, x2, y2, text))
        return result
    
    # Sort labels by angle for easier processing
    sorted_labels = sorted(external_labels, key=lambda x: x[0])
    positioned = []
    
    for i, (angle, x1, y1, text) in enumerate(sorted_labels):
        # Start with default position
        radius = 1.4
        x2 = radius * np.cos(np.radians(angle))
        y2 = radius * np.sin(np.radians(angle))
        
        # Check for overlaps with already positioned labels
        overlap_found = True
        attempts = 0
        max_attempts = 10
        
        while overlap_found and attempts < max_attempts:
            overlap_found = False
            
            for _, _, _, prev_x2, prev_y2, _ in positioned:
                distance = np.sqrt((x2 - prev_x2)**2 + (y2 - prev_y2)**2)
                if distance < min_distance:
                    overlap_found = True
                    break
            
            if overlap_found:
                # Try increasing radius
                radius = min(radius + 0.15, max_radius)
                x2 = radius * np.cos(np.radians(angle))
                y2 = radius * np.sin(np.radians(angle))
                
                # If still overlapping, try small angular adjustments
                if attempts > 5:
                    angle_offset = 15 * (1 if attempts % 2 == 0 else -1)
                    adjusted_angle = angle + angle_offset
                    x2 = radius * np.cos(np.radians(adjusted_angle))
                    y2 = radius * np.sin(np.radians(adjusted_angle))
            
            attempts += 1
        
        positioned.append((angle, x1, y1, x2, y2, text))
    
    return positioned

def load_genome_data(genomesize_file):
    """Load genome size data from file."""
    print(f"Reading genome size file: {genomesize_file}")
    
    try:
        genome_data = pd.read_csv(
            genomesize_file, 
            sep='\t', 
            names=['ID', 'genome_size', 'mutation_rate', 'species_id']
        )
        print(f"Loaded data for {len(genome_data)} samples")
        return genome_data
    except Exception as e:
        print(f"Error reading genome size file: {e}")
        sys.exit(1)

def process_bed_file(bed_file, genome_size, min_size_bp, max_divergence):
    """Process a single BED file and calculate class proportions."""
    print(f"Processing {bed_file}")
    
    # Column names for BED file
    bed_columns = ['Scaffold', 'Start', 'End', 'TE', 'Hit_size', 
                   'Orientation', 'Class', 'Family', 'Divergence', 'RM_ID']
    
    try:
        # Read BED file
        bed_data = pd.read_csv(bed_file, sep='\t', names=bed_columns)
        print(f"  Loaded {len(bed_data)} repeat elements")
        
        # Filter by divergence
        if max_divergence > 0:
            bed_data = bed_data[bed_data['Divergence'] <= max_divergence]
            print(f"  After divergence filter (<= {max_divergence}): {len(bed_data)} elements")
        
        # Filter by minimum size
        if not min_size_bp:  # If flag not set, exclude < 100bp
            bed_data = bed_data[bed_data['Hit_size'] >= 100]
            print(f"  After size filter (>= 100bp): {len(bed_data)} elements")
        
        # Calculate proportions for each class
        class_proportions = bed_data.groupby('Class')['Hit_size'].sum() / genome_size
        
        return class_proportions.to_dict()
        
    except FileNotFoundError:
        print(f"  Warning: File {bed_file} not found. Skipping.")
        return {}
    except Exception as e:
        print(f"  Error processing {bed_file}: {e}")
        return {}

def create_pie_chart(proportions_dict, sample_id, output_prefix):
    """Create pie chart for a single sample."""
    
    # Define class order and colors
    class_order = ['DNA', 'RC', 'LINE', 'SINE', 'LTR', 'Unknown', 'Satellite', 'Simple_repeat']
    class_colors = {
        'DNA': '#FF0000',           # Red
        'RC': '#800080',            # Purple
        'LINE': '#0000FF',          # Blue
        'SINE': '#87CEEB',          # Light blue
        'LTR': '#008000',           # Green
        'Unknown': '#2F2F2F',       # Dark grey
        'Satellite': '#808080',     # Lighter grey
        'Simple_repeat': '#D3D3D3', # Lightest grey
        'Unmasked': '#000000'       # Black
    }
    
    # Prepare data in correct order
    ordered_proportions = []
    ordered_labels = []
    ordered_colors = []
    
    total_masked = 0
    
    # Add classes in specified order
    for class_name in class_order:
        if class_name in proportions_dict:
            prop = proportions_dict[class_name]
            ordered_proportions.append(prop)
            ordered_labels.append(class_name)
            ordered_colors.append(class_colors[class_name])
            total_masked += prop
    
    # Add unmasked proportion
    unmasked_prop = max(0, 1.0 - total_masked)
    ordered_proportions.append(unmasked_prop)
    ordered_labels.append('Unmasked')
    ordered_colors.append(class_colors['Unmasked'])
    
    # Create figure with extra space for external labels
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Create pie chart without labels initially
    wedges, texts = ax.pie(
        ordered_proportions,
        labels=None,
        colors=ordered_colors,
        autopct=None,  # We'll handle percentages manually
        startangle=90,  # Start at 12 o'clock
        textprops={'fontsize': 12}
    )
    
    # Collect external labels first to handle overlaps
    external_labels = []
    internal_labels = []
    
    for i, (wedge, label, prop) in enumerate(zip(wedges, ordered_labels, ordered_proportions)):
        percentage = prop * 100
        angle = (wedge.theta1 + wedge.theta2) / 2
        
        if percentage >= 10:
            # Large wedges: label inside
            x = 0.6 * np.cos(np.radians(angle))
            y = 0.6 * np.sin(np.radians(angle))
            font_color = 'white' if label in ['Unmasked', 'Unknown'] else 'black'
            internal_labels.append((x, y, f'{label}\n{percentage:.1f}%', font_color))
        else:
            # Small wedges: collect for external positioning
            x1 = 1.1 * np.cos(np.radians(angle))
            y1 = 1.1 * np.sin(np.radians(angle))
            external_labels.append((angle, x1, y1, f'{label}\n{percentage:.1f}%'))
    
    # Add internal labels
    for x, y, text, color in internal_labels:
        ax.text(x, y, text, ha='center', va='center', 
               fontsize=10, fontweight='bold', color=color)
    
    # Position external labels with overlap avoidance
    if external_labels:
        positioned_labels = resolve_label_overlaps(external_labels)
        
        for angle, x1, y1, x2, y2, text in positioned_labels:
            # Draw line from wedge to label
            ax.plot([x1, x2], [y1, y2], 'k-', linewidth=1)
            
            # Add label with percentage
            ax.text(x2, y2, text, ha='center', va='center', 
                   fontsize=9, bbox=dict(boxstyle='round,pad=0.3', 
                   facecolor='white', alpha=0.9, edgecolor='gray'))
    
    ax.set_title(f'{sample_id}', fontsize=16, fontweight='bold', pad=20)
    ax.axis('equal')
    
    # Add legend for individual chart
    legend_elements = [plt.Rectangle((0,0),1,1, facecolor=color) 
                      for label, color in zip(ordered_labels, ordered_colors)]
    ax.legend(legend_elements, ordered_labels, loc='center left', 
             bbox_to_anchor=(1.05, 0.5), fontsize=10)
    
    # Adjust layout to accommodate external labels and legend
    plt.tight_layout()
    plt.subplots_adjust(left=0.1, right=0.75, top=0.9, bottom=0.1)
    
    # Save individual pie chart
    plt.savefig(f'{output_prefix}_{sample_id}_piechart.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_{sample_id}_piechart.svg', bbox_inches='tight')
    plt.close()
    
    return ordered_labels, ordered_proportions, ordered_colors

def create_combined_pie_charts(all_results, output_prefix):
    """Create a single multipanel figure with all pie charts."""
    
    n_samples = len(all_results)
    
    # Calculate subplot layout
    cols = min(3, n_samples)
    rows = int(np.ceil(n_samples / cols))
    
    # Create figure with extra space for legend
    fig_width = cols * 6 + 4  # Extra space for legend
    fig_height = rows * 6
    fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(fig_width, fig_height))
    
    # Handle single subplot case
    if n_samples == 1:
        axes = [axes]
    elif rows == 1:
        axes = axes.flatten()
    else:
        axes = axes.flatten()
    
    # Define class order and colors (same as individual charts)
    class_order = ['DNA', 'RC', 'LINE', 'SINE', 'LTR', 'Unknown', 'Satellite', 'Simple_repeat', 'Unmasked']
    class_colors = {
        'DNA': '#FF0000', 'RC': '#800080', 'LINE': '#0000FF', 'SINE': '#87CEEB',
        'LTR': '#008000', 'Unknown': '#2F2F2F', 'Satellite': '#808080', 
        'Simple_repeat': '#D3D3D3', 'Unmasked': '#000000'
    }
    
    legend_elements = []
    legend_labels = []
    
    # Create pie charts for each sample
    for i, (sample_id, (labels, proportions, colors)) in enumerate(all_results.items()):
        ax = axes[i]
        
        # Create pie chart
        wedges, texts = ax.pie(
            proportions,
            colors=colors,
            autopct=None,  # Handle percentages manually
            startangle=90,
            textprops={'fontsize': 8}
        )
        
        # Handle labels for multipanel figure
        external_labels = []
        internal_labels = []
        
        for j, (wedge, label, prop) in enumerate(zip(wedges, labels, proportions)):
            percentage = prop * 100
            angle = (wedge.theta1 + wedge.theta2) / 2
            
            if percentage >= 10:
                # Large wedges: label inside
                x = 0.6 * np.cos(np.radians(angle))
                y = 0.6 * np.sin(np.radians(angle))
                font_color = 'white' if label in ['Unmasked', 'Unknown'] else 'black'
                internal_labels.append((x, y, f'{label}\n{percentage:.1f}%', font_color))
            else:
                # Small wedges: collect for external positioning
                x1 = 1.05 * np.cos(np.radians(angle))
                y1 = 1.05 * np.sin(np.radians(angle))
                external_labels.append((angle, x1, y1, f'{label}\n{percentage:.1f}%'))
        
        # Add internal labels
        for x, y, text, color in internal_labels:
            ax.text(x, y, text, ha='center', va='center', 
                   fontsize=8, fontweight='bold', color=color)
        
        # Position external labels with overlap avoidance (smaller spacing for multipanel)
        if external_labels:
            positioned_labels = resolve_label_overlaps(external_labels, min_distance=0.25, max_radius=1.4)
            
            for angle, x1, y1, x2, y2, text in positioned_labels:
                # Draw line from wedge to label
                ax.plot([x1, x2], [y1, y2], 'k-', linewidth=0.8)
                
                # Add label with percentage
                ax.text(x2, y2, text, ha='center', va='center', 
                       fontsize=7, bbox=dict(boxstyle='round,pad=0.2', 
                       facecolor='white', alpha=0.9, edgecolor='gray'))
        
        ax.set_title(f'{sample_id}', fontsize=12, fontweight='bold')
        ax.axis('equal')
        
        # Collect legend elements from first chart
        if i == 0:
            for label, color in zip(labels, colors):
                legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=color))
                legend_labels.append(label)
    
    # Hide unused subplots
    for j in range(n_samples, len(axes)):
        axes[j].axis('off')
    
    # Add single legend to the right of the figure
    fig.legend(legend_elements, legend_labels, loc='center right', 
              bbox_to_anchor=(0.95, 0.5), fontsize=12)
    
    plt.tight_layout()
    plt.subplots_adjust(right=0.8)  # Make room for legend
    
    # Save combined figure
    plt.savefig(f'{output_prefix}_all_samples_piecharts.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_all_samples_piecharts.svg', bbox_inches='tight')
    plt.close()
    
    print(f"Combined pie chart saved as {output_prefix}_all_samples_piecharts.png/svg")

def save_proportion_tables(all_results, output_prefix):
    """Save proportion tables for all samples."""
    
    # Create summary table
    summary_data = []
    
    for sample_id, (labels, proportions, colors) in all_results.items():
        row = {'Sample_ID': sample_id}
        for label, prop in zip(labels, proportions):
            row[label] = f"{prop * 100:.2f}%"
        summary_data.append(row)
    
    # Convert to DataFrame and save
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(f'{output_prefix}_proportion_summary.csv', index=False)
    
    print(f"Proportion summary saved as {output_prefix}_proportion_summary.csv")
    
    # Save individual tables
    for sample_id, (labels, proportions, colors) in all_results.items():
        individual_df = pd.DataFrame({
            'Class': labels,
            'Proportion': proportions,
            'Percentage': [f"{p * 100:.2f}%" for p in proportions]
        })
        individual_df.to_csv(f'{output_prefix}_{sample_id}_proportions.csv', index=False)

def main():
    """Main function."""
    args = get_args()
    
    # Load genome data
    genome_data = load_genome_data(args.genomesize)
    
    # Process each sample
    all_results = {}
    output_prefix = "repeat_analysis"
    
    print(f"\nProcessing {len(genome_data)} samples...")
    print(f"Divergence threshold: {args.divergence}")
    print(f"Include hits < 100bp: {args.minimum100bp}")
    
    for _, row in genome_data.iterrows():
        sample_id = row['ID']
        genome_size = row['genome_size']
        
        # Construct BED file name
        bed_file = f"{sample_id}_rm.bed"
        
        # Process BED file
        proportions = process_bed_file(
            bed_file, 
            genome_size, 
            args.minimum100bp, 
            args.divergence
        )
        
        if proportions:  # Only process if we got data
            # Create individual pie chart
            labels, props, colors = create_pie_chart(proportions, sample_id, output_prefix)
            all_results[sample_id] = (labels, props, colors)
            print(f"  Individual pie chart saved for {sample_id}")
    
    if all_results:
        # Create combined pie charts
        print(f"\nCreating combined pie chart with {len(all_results)} samples...")
        create_combined_pie_charts(all_results, output_prefix)
        
        # Save proportion tables
        print("Saving proportion tables...")
        save_proportion_tables(all_results, output_prefix)
        
        print(f"\nAnalysis complete! Processed {len(all_results)} samples.")
    else:
        print("No data was successfully processed. Please check your input files.")

if __name__ == "__main__":
    main()
