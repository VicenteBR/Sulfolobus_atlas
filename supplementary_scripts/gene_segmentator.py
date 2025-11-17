import os
import pandas as pd
import pybedtools
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator

def gff_to_bed(gff_file):
    """Convert a GFF file to a BED format for easier intersection."""
    bed_data = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] == "gene":  # Extract only gene entries
                chrom, start, end, strand, attributes = fields[0], int(fields[3]) - 1, int(fields[4]), fields[6], fields[8]
                gene_id = [x for x in attributes.split(";") if x.startswith("Name=")]
                if gene_id:
                    gene_id = gene_id[0].replace("Name=", "")
                    bed_data.append([chrom, start, end, gene_id, ".", strand])
    return pybedtools.BedTool(bed_data)

def feature_to_bed(feature_gff, feature_name):
    """Convert a feature GFF file to BED format."""
    bed_data = []
    with open(feature_gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom, start, end, strand = fields[0], int(fields[3]) - 1, int(fields[4]), fields[6]
            bed_data.append([chrom, start, end, feature_name, ".", strand])
    return pybedtools.BedTool(bed_data)

def classify_feature_position(feature_start, feature_end, gene_start, gene_end, strand):
    """Classify the feature position within the gene, considering strand orientation."""
    if strand == "-":
        gene_start, gene_end = gene_end, gene_start
        feature_start, feature_end = feature_end, feature_start

    gene_length = gene_end - gene_start
    feature_mid = (feature_start + feature_end) / 2
    relative_position = (feature_mid - gene_start) / gene_length

    if relative_position < 0.25:
        return "5'" if strand == "+" else "3'"
    elif 0.25 <= relative_position < 0.75:
        return "internal"
    else:
        return "3" if strand == "+" else "5"

def plot_gene_with_features(gene_start, gene_end, features, gene_id, feature_names, feature_colors, strand):
    """Plot a gene with features and save as an image, using a bar to indicate gene direction."""
    plt.figure(figsize=(12, 6))  # Increased height for better spacing

    # Plot the gene as a bar
    if strand == "+":
        plt.plot([gene_start, gene_end], [0, 0], color="black", lw=15)  # Main gene bar
    else:
        plt.plot([gene_end, gene_start], [0, 0], color="black", lw=15)  # Main gene bar

    # Plot features
    feature_positions = {}  # Dictionary to track the y-position for each feature
    y_offset = 1
    for feature in features:
        feature_start, feature_end, feature_name = feature
        if feature_name not in feature_positions:
            feature_positions[feature_name] = y_offset
            y_offset += 1

        feature_length = feature_end - feature_start

        # Plot small features as vertical dashed lines if their length is less than 15
        if feature_length < 15:
            plt.plot([feature_start, feature_start], 
                     [feature_positions[feature_name] - 0.4, feature_positions[feature_name] + 0.4], 
                     color=feature_colors[feature_name], ls="--", lw=2)
        else:
            plt.plot([feature_start, feature_end], 
                     [feature_positions[feature_name], feature_positions[feature_name]], 
                     color=feature_colors[feature_name], lw=8)

    # Add (+) or (-) to gene ID in the plot title
    strand_indicator = "+" if strand == "+" else "-"
    plt.title(f"Gene: {gene_id} ({strand_indicator})")
    plt.yticks([])  # Hide y-axis
    plt.xlabel("Genomic Position")

    # Set x-axis limits to the start and end of the gene (no padding)
    plt.xlim(min(gene_start, gene_end), max(gene_start, gene_end))

    # Ensure the x-axis is not in scientific notation and does not display negative values
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True, prune='lower'))  # Prune lower side to avoid negative values
    plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{int(x)}'))  # Ensure no scientific notation

    # Add legend with colors for each feature
    handles = [mpatches.Patch(color=feature_colors[name], label=name) for name in feature_names]
    plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left')

    # Save plot to 'plots' folder
    plot_dir = "plots"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/{gene_id}_{strand_indicator}_features.png")
    plt.close()

def create_feature_matrices(main_gff, feature_folder):
    """Create both a binary presence-absence matrix and a feature count table, considering strand orientation."""
    gene_bed = gff_to_bed(main_gff)
    
    # Get all GFF files in the features folder
    feature_gff_files = [os.path.join(feature_folder, f) for f in os.listdir(feature_folder) if f.endswith(".gff")]
    feature_names = [os.path.basename(f).replace(".gff", "") for f in feature_gff_files]
    
    genes = {entry[3]: [0] * len(feature_gff_files) for entry in gene_bed}
    gene_feature_counts = {entry[3]: [0] * len(feature_gff_files) for entry in gene_bed}
    gene_feature_positions = {entry[3]: [[] for _ in range(len(feature_gff_files))] for entry in gene_bed}
    gene_features = {entry[3]: [] for entry in gene_bed}

    feature_colors = {name: plt.cm.tab10(idx % 10) for idx, name in enumerate(feature_names)}

    for idx, feature_file in enumerate(feature_gff_files):
        feature_bed = feature_to_bed(feature_file, feature_names[idx])
        
        # Find intersections, considering strand
        intersections = gene_bed.intersect(feature_bed, wa=True, wb=True, s=True)  # `-s` ensures strand match

        for entry in intersections:
            gene_id = entry[3]  # Gene ID is in the 4th column of gene BED
            feature_start = int(entry[7])  # Corrected: Feature start in the 8th column of feature BED
            feature_end = int(entry[8])  # Corrected: Feature end in the 9th column of feature BED
            gene_start = int(entry[1])  # Gene start in the 2nd column of gene BED
            gene_end = int(entry[2])  # Gene end in the 3rd column of gene BED
            gene_strand = entry[5]  # The strand information of the gene

            if gene_id in genes:
                genes[gene_id][idx] = 1  # Binary presence
                gene_feature_counts[gene_id][idx] += 1  # Count occurrences
                position = classify_feature_position(feature_start, feature_end, gene_start, gene_end, gene_strand)
                gene_feature_positions[gene_id][idx].append(position)
                gene_features[gene_id].append((feature_start, feature_end, feature_names[idx]))

    # Plot each gene
    for gene_id, features in gene_features.items():
        # Find the gene start and end by searching in the gene_bed
        for entry in gene_bed:
            if entry[3] == gene_id:
                gene_start = int(entry[1])
                gene_end = int(entry[2])
                gene_strand = entry[5]  # Strand information
                break
        plot_gene_with_features(gene_start, gene_end, features, gene_id, feature_names, feature_colors, gene_strand)

    # Create DataFrames
    binary_matrix_df = pd.DataFrame.from_dict(genes, orient="index", columns=feature_names)
    binary_matrix_df.index.name = "Gene"
    
    feature_counts_df = pd.DataFrame.from_dict(gene_feature_counts, orient="index", columns=feature_names)
    feature_counts_df.index.name = "Gene"
    
    feature_positions_df = pd.DataFrame.from_dict(gene_feature_positions, orient="index", columns=feature_names)
    feature_positions_df.index.name = "Gene"
    
    return binary_matrix_df, feature_counts_df, feature_positions_df

if __name__ == "__main__":
    # Get user input for the gene annotation file
    main_gff = input("Enter the path to the gene annotation GFF file: ").strip()
    feature_folder = "features"  # Folder containing feature GFF files

    # Check if feature folder exists
    if not os.path.isdir(feature_folder):
        print(f"Error: Feature folder '{feature_folder}' not found.")
        exit(1)

    binary_matrix, feature_counts, feature_positions = create_feature_matrices(main_gff, feature_folder)
    
    binary_matrix.to_csv("feature_matrix.csv")
    feature_counts.to_csv("feature_counts.csv")
    feature_positions.to_csv("feature_positions.csv")

    print("Feature matrix saved to 'feature_matrix.csv'")
    print("Feature counts saved to 'feature_counts.csv'")
    print("Feature positions saved to 'feature_positions.csv'")
    print(binary_matrix)
