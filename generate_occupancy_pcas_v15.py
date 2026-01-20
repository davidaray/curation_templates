import os
import argparse
import pandas as pd
from collections import defaultdict
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import multiprocessing as mp
import gzip

# Use 'Agg' backend for headless environments
matplotlib.use('Agg')

# List of TE classes and families, including Satellite
te_classes = ["LINE", "SINE", "LTR", "DNA", "RC", "Satellite"]

# Function to load mapping file for species and their taxonomic families
def load_species_family_mapping(mapping_file):
    species_family_map = {}
    with open(mapping_file, 'r') as file:
        next(file)  # Skip header
        for line in file:
            line = line.strip()
            if not line:
                continue
            columns = line.split('\t')
            if len(columns) < 6:
                print(f"Warning: Skipping malformed line in mapping file: {line}")
                continue
            tax_family = columns[0]
            species_id = columns[2]
            species_family_map[species_id] = tax_family
    return species_family_map

# Function to sanitize filenames by replacing spaces with underscores
def sanitize_filename(filename):
    return filename.replace(' ', '_')

# Function to extract class and family data from .out.gz files
def extract_class_family_data(file_path, species_family_map):
    te_class_data = defaultdict(int)
    te_family_data = {cls: defaultdict(int) for cls in te_classes}
    species = os.path.basename(file_path).split('.out.gz')[0]
    tax_family = species_family_map.get(species, "Unknown")

    with gzip.open(file_path, 'rt') as file:
        for _ in range(3):
            next(file)
        for line in file:
            columns = line.strip().split()
            if len(columns) > 10:
                te_class_family = columns[10]
                te_class = te_class_family.split('/')[0]
                te_family = te_class_family.split('/')[1] if '/' in te_class_family else "Unknown"
                try:
                    start = int(columns[5])
                    end = int(columns[6])
                    occupancy = end - start + 1
                    if te_class in te_classes:
                        te_class_data[te_class] += occupancy
                        te_family_data[te_class][te_family] += occupancy
                except ValueError:
                    print(f"Skipping malformed line in {file_path}: {line}")
    return species, tax_family, te_class_data, te_family_data

# Function to build matrices using multiprocessing
def build_matrices(file_list, species_family_map, num_proc):
    with mp.Pool(num_proc) as pool:
        results = pool.starmap(
            extract_class_family_data,
            [(file_path, species_family_map) for file_path in file_list]
        )
    class_matrix = pd.DataFrame()
    family_matrices = {cls: pd.DataFrame() for cls in te_classes}
    for species, tax_family, te_class_data, te_family_data in results:
        class_matrix = pd.concat([class_matrix, pd.DataFrame(te_class_data, index=[species])], axis=0).fillna(0)
        for te_class, family_occupancys in te_family_data.items():
            family_matrices[te_class] = pd.concat([family_matrices[te_class], pd.DataFrame(family_occupancys, index=[species])], axis=0).fillna(0)
    return class_matrix, family_matrices

# Function to save matrices
def save_matrices(class_matrix, family_matrices, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    class_matrix.to_csv(sanitize_filename(os.path.join(output_dir, 'class_PCA_occupancy_matrix.tsv')), sep='\t')
    for te_class, df in family_matrices.items():
        df.to_csv(sanitize_filename(os.path.join(output_dir, f'{te_class}_family_PCA_occupancy_matrix.tsv')), sep='\t')

def run_pca_and_plot(matrix, output_prefix, species_family_map):
    custom_palette = {
        "Pteropodidae": "#484871",
        "Rhinopomatidae": "#BF9000",
        "Megadermatidae": "#FFF2CC",
        "Craseonycteridae": "#806000",
        "Rhinonycteridae": "#FFD966",
        "Hipposideridae": "#CE9B50",
        "Rhinolophidae": "#FFC000",
        "Emballonuridae": "#9B9BD4",
        "Nycteridae": "#6262BC",
        "Mystacinidae": "#FF8877",
        "Thyropteridae": "#F35A2B",
        "Furipteridae": "#F6B26B",
        "Noctilionidae": "#F8CBAD",
        "Mormoopidae": "#F87523",
        "Phyllostomidae": "#F9CB9C",
        "Myzopodidae": "#D7F4EE",
        "Natalidae": "#267D73",
        "Molossidae": "#6FC3B6",
        "Miniopteridae": "#2C6354",
        "Vespertilionidae": "#D9EAD3",
        "Cistugidae": "#38A884",
        "Unknown": "#CCCCCC"  # Default color for unknown taxonomic families
    }
    # Define custom markers for specific families
    custom_markers = {
        "Pteropodidae": "o",      # Circle
        "Rhinopomatidae": "s",    # Square
        "Megadermatidae": "D",    # Diamond
        "Craseonycteridae": "^",  # Upward triangle
        "Rhinonycteridae": "v",   # Downward triangle
        "Hipposideridae": "<",    # Left triangle
        "Rhinolophidae": ">",     # Right triangle
        "Emballonuridae": "p",    # Pentagon
        "Nycteridae": "H",        # Hexagon
        "Mystacinidae": "*",      # Star
        "Thyropteridae": "X",     # X (filled)
        "Furipteridae": "8",      # Octagon
        "Noctilionidae": "P",     # Plus (filled)
        "Mormoopidae": "d",       # Thin diamond
        "Phyllostomidae": "h",    # Hexagon alternative
        "Myzopodidae": "o",       # Circle
        "Natalidae": "s",         # Square
        "Molossidae": "d",        # Thin diamond
        "Miniopteridae": "p",     # Pentagon
        "Vespertilionidae": "H",  # Hexagon
        "Cistugidae": "X",        # X (filled)
        "Unknown": "o"            # Default circle
    }

    # Create nicer label for titles
    base_name = os.path.basename(output_prefix)
    if base_name == "class_matrix":
        label = "TE Class Insertion Occupancy"
    elif base_name.endswith("_family_matrix"):
        class_name = base_name.replace("_family_matrix", "")
        label = f"{class_name} Insertion Occupancy"
    else:
        label = base_name

    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(matrix)
    
    # Create DataFrame for PCA results
    pca_df = pd.DataFrame(pca_result, index=matrix.index, columns=['PC1', 'PC2'])
    pca_df['Taxonomic_Family'] = pca_df.index.map(lambda x: species_family_map.get(x, 'Unknown'))

    # Save PCA results
    pca_df.to_csv(f'{output_prefix}_PCA_occupancy_scores.tsv', sep='\t')
    
    # Scatter plot
    plt.figure(figsize=(10, 7))
    scatter_plot = sns.scatterplot(
        data=pca_df,
        x='PC1', y='PC2', hue='Taxonomic_Family', style='Taxonomic_Family',
        palette=custom_palette, markers=custom_markers, s=100, alpha=0.7, edgecolor="black"
    )
    plt.title(f'PCA Plot for {label}')
    plt.xlabel(f'PC1 - {pca.explained_variance_ratio_[0]*100:.2f}%')
    plt.ylabel(f'PC2 - {pca.explained_variance_ratio_[1]*100:.2f}%')

    # Customize legend
    handles, labels = scatter_plot.get_legend_handles_labels()
    sorted_labels, sorted_handles = zip(*sorted(zip(labels, handles)))
    legend = plt.legend(
        handles=sorted_handles, labels=sorted_labels,
        loc='center left', bbox_to_anchor=(1, 0.5),
        title='Taxonomic Family', frameon=True,
        handlelength=4.0
    )
    marker_size = 100
    for handle in legend.legend_handles:
        handle.set_sizes([marker_size])
        handle.set_edgecolor('black')

    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(f'{output_prefix}_PCA_occupancy_plot.png', bbox_inches='tight')
    plt.close()
    
    # Biplot with arrows
    plt.figure(figsize=(10, 7))
    sns.scatterplot(
        data=pca_df,
        x='PC1', y='PC2', hue='Taxonomic_Family', style='Taxonomic_Family',
        palette=custom_palette, markers=custom_markers, s=100, alpha=0.7, edgecolor="black"
    )
    for i, (pc1, pc2) in enumerate(zip(pca.components_[0], pca.components_[1])):
        plt.arrow(
            0, 0, pc1 * max(pca_df['PC1']), pc2 * max(pca_df['PC2']),
            color='black', alpha=0.5, head_width=0.02, head_length=0.02
        )
        plt.text(
            pc1 * max(pca_df['PC1']) * 1.1, pc2 * max(pca_df['PC2']) * 1.1,
            matrix.columns[i], color='black', ha='center', va='center'
        )
    plt.title(f'PCA Biplot for {label}')
    plt.xlabel(f'PC1 - {pca.explained_variance_ratio_[0]*100:.2f}%')
    plt.ylabel(f'PC2 - {pca.explained_variance_ratio_[1]*100:.2f}%')
    
    # Customize legend
    handles, labels = scatter_plot.get_legend_handles_labels()
    sorted_labels, sorted_handles = zip(*sorted(zip(labels, handles)))
    legend = plt.legend(
        handles=sorted_handles, labels=sorted_labels,
        loc='center left', bbox_to_anchor=(1, 0.5),
        title='Taxonomic Family', frameon=True,
        handlelength=4.0
    )
    for handle in legend.legend_handles:
        handle.set_sizes([marker_size])
        handle.set_edgecolor('black')

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_PCA_occupancy_biplot.png')
    plt.close()

    # Save PCA results to file
    pca_df.drop(columns='Taxonomic_Family').to_csv(f'{output_prefix}_PCA_occupancy_result.tsv', sep='\t')
    pca_df[['PC1', 'PC2', 'Taxonomic_Family']].to_csv(f'{output_prefix}_occupancy_coordinates.tsv', sep='\t')

    # Analyze and save PCA loadings
    loadings = pca.components_
    loadings_df = pd.DataFrame(loadings.T, index=matrix.columns, columns=['PC1', 'PC2'])
    loadings_df.to_csv(f'{output_prefix}_PCA_occupancy_loadings.tsv', sep='\t')

    print(f'PCA Loadings for {output_prefix} saved to {output_prefix}_PCA_occupancy_loadings.tsv')

# Main function
def main():
    parser = argparse.ArgumentParser(description='Generate PCA matrices and plots from RepeatMasker .out.gz files.')
    parser.add_argument('--out_dir', required=True, help='Directory containing .out.gz files.')
    parser.add_argument('--mapping_file', required=True, help='File mapping species to taxonomic families.')
    parser.add_argument('--output_dir', required=True, help='Directory to save output matrices and PCA plots.')
    parser.add_argument('--num_proc', type=int, default=1, help='Number of processors to use for multiprocessing.')
    parser.add_argument('--skip_calculations', action='store_true', help='Skip matrix calculations and go straight to plotting.')
    
    args = parser.parse_args()
    
    # Load species family mapping
    species_family_map = load_species_family_mapping(args.mapping_file)
    
    # Get list of .out.gz files
    file_list = [os.path.join(args.out_dir, f) for f in os.listdir(args.out_dir) if f.endswith('.out.gz')]
    
    # Check if matrices already exist
    class_matrix_path = os.path.join(args.output_dir, 'class_PCA_occupancy_matrix.tsv')
    family_matrix_paths = {te_class: os.path.join(args.output_dir, f'{te_class}_family_PCA_occupancy_matrix.tsv') for te_class in te_classes}
    
    matrices_exist = os.path.exists(class_matrix_path) and all(os.path.exists(family_matrix_paths[te_class]) for te_class in te_classes)
    
    if not args.skip_calculations and not matrices_exist:
        # Build matrices using multiprocessing if they don't exist
        class_matrix, family_matrices = build_matrices(file_list, species_family_map, args.num_proc)
        
        # Save matrices
        save_matrices(class_matrix, family_matrices, args.output_dir)
    else:
        # Load existing matrices
        class_matrix = pd.read_csv(class_matrix_path, sep='\t', index_col=0)
        family_matrices = {
            te_class: pd.read_csv(family_matrix_paths[te_class], sep='\t', index_col=0) for te_class in te_classes
        }
    
    # Run PCA for class matrix
    run_pca_and_plot(class_matrix, os.path.join(args.output_dir, 'class_matrix'), species_family_map)
    
    # Run PCA for each family matrix
    for te_class, family_matrix in family_matrices.items():
        if not family_matrix.empty:
            run_pca_and_plot(family_matrix, os.path.join(args.output_dir, f'{te_class}_family_matrix'), species_family_map)

if __name__ == '__main__':
    main()

