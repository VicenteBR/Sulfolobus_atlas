import os
import pybedtools
from Bio import SeqIO

def gff_to_bed(gff_file):
    """Convert a GFF file to a BED format for easier intersection."""
    bed_data = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] == "gene":
                chrom, start, end, strand, attributes = (
                    fields[0], int(fields[3]) - 1, int(fields[4]), fields[6], fields[8]
                )
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
        return "5" if strand == "+" else "3"
    elif 0.25 <= relative_position < 0.75:
        return "internal"
    else:
        return "3" if strand == "+" else "5"

def calculate_gc_content(sequence):
    """Calculate GC content of a given sequence."""
    gc_count = sum(1 for base in sequence if base in "GCgc")
    return round((gc_count / len(sequence)) * 100, 2) if sequence else 0.0

def extract_feature_sequences(fasta_file, bed_entries):
    """Extract sequences from the genome FASTA based on feature locations and separate them by classification."""
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    fasta_output = {"3": [], "5": [], "internal": []}
    
    for entry in bed_entries:
        chrom, start, end, feature_name, gene_id, strand, position = entry
        start, end = int(start), int(end)

        if chrom in genome:
            sequence = genome[chrom].seq[start:end]
            if strand == "-":
                sequence = sequence.reverse_complement()

            gc_content = calculate_gc_content(sequence)
            fasta_output[position].append(f">{gene_id}|{feature_name}|{chrom}:{start}-{end}({strand})|GC:{gc_content}%\n{sequence}\n")
    
    output_dir = "fasta_outputs"
    os.makedirs(output_dir, exist_ok=True)
    for category, sequences in fasta_output.items():
        if sequences:
            file_path = os.path.join(output_dir, f"{category}_features.fasta")
            with open(file_path, "w") as f:
                f.writelines(sequences)
            print(f"Saved: {file_path}")

def create_feature_matrices(main_gff, feature_folder, fasta_file):
    """Extract features and classify them based on their position."""
    gene_bed = gff_to_bed(main_gff)
    feature_gff_files = [os.path.join(feature_folder, f) for f in os.listdir(feature_folder) if f.endswith(".gff")]
    feature_names = [os.path.basename(f).replace(".gff", "") for f in feature_gff_files]
    gene_features = {entry[3]: [] for entry in gene_bed}
    
    for idx, feature_file in enumerate(feature_gff_files):
        feature_bed = feature_to_bed(feature_file, feature_names[idx])
        intersections = gene_bed.intersect(feature_bed, wa=True, wb=True, s=True)
        
        for entry in intersections:
            gene_id = entry[3]
            feature_start = int(entry[7])
            feature_end = int(entry[8])
            gene_start = int(entry[1])
            gene_end = int(entry[2])
            gene_strand = entry[5]

            if gene_id in gene_features:
                position = classify_feature_position(feature_start, feature_end, gene_start, gene_end, gene_strand)
                gene_features[gene_id].append((entry[0], feature_start, feature_end, feature_names[idx], gene_id, gene_strand, position))
    
    if fasta_file:
        extract_feature_sequences(fasta_file, [entry for values in gene_features.values() for entry in values])

if __name__ == "__main__":
    main_gff = input("Enter the path to the gene annotation GFF file: ").strip()
    fasta_file = input("Enter the path to the genome FASTA file: ").strip()
    feature_folder = "features"
    
    if not os.path.isdir(feature_folder):
        print(f"Error: Feature folder '{feature_folder}' not found.")
        exit(1)

    create_feature_matrices(main_gff, feature_folder, fasta_file)