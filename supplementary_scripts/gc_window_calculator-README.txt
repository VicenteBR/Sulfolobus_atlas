README — gc_window_calculator.py

Overview
This script calculates the GC content of DNA sequences using sliding windows. It reads a FASTA file, scans each sequence with a 50 bp window moving in 25 bp steps, computes GC percentage for every window, and outputs the results to a tab-delimited file.

Features
- Reads FASTA files via BioPython.
- Sliding window calculation:
  - Window size: 50 bp
  - Step size: 25 bp
- Computes GC content for each window.
- Outputs a table containing Chromosome, Start, End, GC_Content.
- Ensures windows do not exceed sequence boundaries.
- Works on sequences of any length.

Usage
python gc_window_calculator.py

You will be prompted for the FASTA file path:
Enter the path to your FASTA file:

Output
The script generates:
gc_content_windows.txt

Each line has the format:
Chromosome    Start    End    GC_Content

Method Summary
1. Parse FASTA records using BioPython.
2. Convert the sequence to uppercase.
3. Iterate over the sequence using:
   - Window size: 50 bp
   - Step: 25 bp
4. Compute GC content:
   GC% = (G_count + C_count) / window_length * 100
5. Save all results in a tab-delimited file.

Dependencies
- Python >= 3.7
- BioPython

Author
José Vicente Gomes Filho
