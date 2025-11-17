README - gc_region_calculator.py

Overview
This script processes gene annotations and genomic feature annotations to classify feature positions along genes, extract their sequences, compute GC content, and generate FASTA outputs grouped into 5′, 3′, and internal categories.

It supports:
- Converting gene and feature GFF files to BED format
- Identifying strand-specific overlaps between genes and arbitrary features
- Classifying each feature relative to the gene (5′, internal, 3′)
- Extracting the corresponding genomic sequence (with reverse-complement for minus strand)
- Calculating GC content for each extracted feature sequence
- Exporting categorized FASTA files

Input Requirements
1. Gene annotation GFF file:
   - Must contain “gene” entries
   - Must include Name= attribute for gene identifiers

2. Genome FASTA file:
   - Must contain sequences matching chromosome identifiers used in GFF files

3. Feature folder ("features/"):
   - Contains multiple .gff files
   - Each file represents one type of genomic feature

Usage
Run the script and provide the required paths when prompted:

Enter the path to the gene annotation GFF file:
Enter the path to the genome FASTA file:

The script expects a directory named “features” containing .gff feature files.

Output Files
A folder named fasta_outputs/ will be created containing:
- 5_features.fasta
- internal_features.fasta
- 3_features.fasta

Each FASTA header includes:
>geneID|featureName|chrom:start-end(strand)|GC:xx.xx%

Method Summary
1. Convert gene GFF → BED (only “gene” entries).
2. Convert each feature GFF → BED (all entries).
3. Perform strand-specific intersection (bedtools -s).
4. For each intersecting feature:
   - Determine midpoint relative to gene
   - Classify as 5′, internal, or 3′ (strand-aware)
5. Extract sequence from genome FASTA (reverse complement if required).
6. Calculate GC%.
7. Write categorized FASTA files.

Dependencies
- Python >= 3.7
- pybedtools
- BioPython
- bedtools installed on system

Author
José Vicente Gomes Filho
