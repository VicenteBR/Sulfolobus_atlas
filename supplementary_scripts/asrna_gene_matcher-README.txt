README — asrna_gene_matcher.py

Overview
This script identifies antisense RNA (asRNA) overlaps with annotated genes using genomic coordinates from GFF files. 
It converts GFF entries to BED format, performs coordinate intersections using pybedtools, and reports which asRNAs overlap each gene on the opposite DNA strand.

Features
- Converts GFF annotations to BED format.
- Supports filtering by feature type (e.g., "gene").
- Extracts Name= attributes to use as feature identifiers.
- Computes overlaps using bedtools via the pybedtools interface.
- Checks opposite-strand orientation to confirm antisense relationships.
- Outputs a clean Gene → asRNAs mapping table.

Usage
python asrna_gene_matcher.py <genes.gff> <asRNAs.gff>

Example:
python asrna_gene_matcher.py annotations/genes.gff annotations/asRNAs.gff

Output:
A tab-separated file:
gene_asrna_matches.tsv

Input Requirements
1. Gene GFF file
   - Must contain feature type "gene" (column 3)
   - Must contain Name= attribute

2. asRNA GFF file
   - May contain only antisense RNAs or mixed feature types
   - Must contain Name= attribute

General steps:
1. Converts GFF to BED
2. Computes intersections with bedtools
3. Filters by opposite strand
4. Outputs Gene → asRNAs table

Dependencies
- Python >= 3.7
- pybedtools
- bedtools
- sys

Author
José Vicente Gomes Filho
