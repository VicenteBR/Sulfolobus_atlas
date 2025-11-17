README — gene_segmentator.py

Overview
This script processes gene annotations and genomic feature annotations to generate:
1. Binary presence–absence matrix (gene × feature)
2. Feature count matrix (number of occurrences per feature per gene)
3. Feature positional matrix (5′ / internal / 3′ classification)
4. Per-gene feature maps, saved as PNG plots

Features
- Converts gene and feature GFFs into BED.
- Matches features to genes on the same strand.
- Classifies features relative to gene (5′, internal, 3′).
- Tracks number and type of features per gene.
- Creates gene-level plots with directional bars and color-coded features.
- Outputs CSV matrices:
  - feature_matrix.csv
  - feature_counts.csv
  - feature_positions.csv

Input Requirements
1. Gene GFF file (must contain gene features with Name= attribute)
2. 'features/' folder containing one or more .gff files (feature type derived from filename)
3. No FASTA file required

Usage
python gene_segmentator.py

You will be prompted for the gene annotation GFF file path.

Outputs
1. Binary matrix: feature_matrix.csv
2. Count matrix: feature_counts.csv
3. Positional matrix: feature_positions.csv
4. Gene feature plots saved in 'plots/' directory.

Method Summary
1. Gene and feature GFF files converted to BED format.
2. Strand-specific intersection used to map features to genes.
3. Features classified by midpoint along gene length.
4. Matrices constructed for presence, counts, and positions.
5. Matplotlib used to visualize each gene with its mapped features.

Dependencies
- Python >= 3.7
- pandas
- matplotlib
- pybedtools
- bedtools installed system-wide

Author
José Vicente Gomes Filho
