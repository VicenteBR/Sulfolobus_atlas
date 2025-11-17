import sys
import pybedtools

def gff_to_bed(gff_file, feature_type=None):
    """
    Convert a GFF file to a BED file.
    If feature_type is provided, only lines where the third column matches it are used.
    Extracts the 'Name=' attribute as the gene (or feature) identifier.
    """
    bed_data = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            # If feature_type is specified, filter by it.
            if feature_type and fields[2] != feature_type:
                continue
            chrom = fields[0]
            start = int(fields[3]) - 1  # BED format uses 0-based start
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]
            # Look for the Name attribute (e.g., Name=gene1)
            feature_id = None
            for attr in attributes.split(";"):
                if attr.startswith("Name="):
                    feature_id = attr.split("=")[1]
                    break
            if feature_id:
                bed_data.append([chrom, start, end, feature_id, ".", strand])
    return pybedtools.BedTool(bed_data)

def main(genes_gff, asrna_gff, output_file="gene_asrna_matches.tsv"):
    # Convert genes and asRNAs to BED.
    # For genes we filter on feature type "gene", for asRNAs we assume the file contains only asRNAs.
    genes_bed = gff_to_bed(genes_gff, feature_type="gene")
    asrna_bed = gff_to_bed(asrna_gff)
    
    # Find intersections.
    # We don't use the strand flag (-s) because we want to verify anti-sense manually.
    intersections = genes_bed.intersect(asrna_bed, wa=True, wb=True)
    
    # Create a dictionary mapping gene IDs to a set of asRNA IDs that overlap them.
    gene_to_asrna = {}
    for entry in intersections:
        # BED file layout:
        # Columns 0-5: gene (chrom, start, end, gene_id, score, strand)
        # Columns 6-11: asRNA (chrom, start, end, asrna_id, score, strand)
        gene_id = entry[3]
        asrna_id = entry[9]
        gene_strand = entry[5]
        asrna_strand = entry[11]
        
        # Check that the asRNA is on the opposite strand to be considered anti-sense.
        if gene_strand != asrna_strand:
            if gene_id not in gene_to_asrna:
                gene_to_asrna[gene_id] = set()
            gene_to_asrna[gene_id].add(asrna_id)
    
    # Save output as a tab-delimited file.
    with open(output_file, "w") as f:
        # Write header.
        f.write("Gene\tasRNAs\n")
        for gene in sorted(gene_to_asrna):
            asrna_list = sorted(gene_to_asrna[gene])
            # Join multiple asRNA IDs with a comma.
            f.write(f"{gene}\t{','.join(asrna_list)}\n")
    
    print(f"Output saved as tab-delimited file: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python match_asrna_genes.py <genes.gff> <asRNAs.gff>")
        sys.exit(1)
    
    genes_gff = sys.argv[1]
    asrna_gff = sys.argv[2]
    main(genes_gff, asrna_gff)
