from Bio import SeqIO

def gc_content(sequence):
    """Calculate GC content percentage in a given DNA sequence."""
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0

def gc_content_windows(fasta_file, window_size=50, step_size=25):
    """Compute GC content using sliding windows from a user-specified FASTA file."""
    results = []

    # Parse the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = record.seq.upper()
        genome_length = len(seq)

        # Sliding window with 25 nt step
        for start in range(0, genome_length - window_size + 1, step_size):
            end = start + window_size
            window_seq = seq[start:end]
            gc_value = gc_content(window_seq)
            results.append((record.id, start, end, gc_value))

    return results

if __name__ == "__main__":
    # Ask user for FASTA file input
    fasta_file = input("Enter the path to your FASTA file: ").strip()

    try:
        # Run the function and save results
        results = gc_content_windows(fasta_file)

        # Output file name
        output_file = "gc_content_windows.txt"

        # Save results to a text file
        with open(output_file, "w") as out_file:
            out_file.write("Chromosome\tStart\tEnd\tGC_Content\n")
            for row in results:
                out_file.write(f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]:.2f}\n")

        print(f"GC content sliding windows saved to '{output_file}'")

    except FileNotFoundError:
        print("Error: File not found. Please check the file path and try again.")
