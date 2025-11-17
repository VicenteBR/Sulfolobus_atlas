README — loci_caller.sh

Overview
This Bash script identifies transcript boundaries by pairing Transcription Start Sites (TSS) with the nearest valid Transcription Termination Site (TTS) within a maximum distance of 400 bp on the same chromosome and strand. 

The script:
1. Loads all TTS entries from tts.gff into an associative array.
2. Reads each TSS entry from tss.gff.
3. For each TSS, searches for the furthest valid TTS within 400 bp in the correct orientation.
4. Outputs transcript intervals as GFF records.

Input Requirements
- tts.gff : GFF file containing TTS entries (type = "TTS")
- tss.gff : GFF file containing TSS entries (type = "TSS")
Both files must follow standard 9-column GFF format.

Output
If a valid TTS is found:
- "+" strand:
    transcript = TSS start → TTS end
- "-" strand:
    transcript = TTS start → TSS start

Output format (tab-delimited):
chr    source    transcript    start    end    .    strand    .    TSS=<attributes>;TTS=<value>

The script prints transcript intervals to STDOUT.
Redirect to file if needed:
./script.sh > transcripts.gff

Dependencies
- Bash 4+ (associative arrays)
- Standard Linux command-line tools
- Optional downstream tools:
  - bedtools
  - awk

Author
José Vicente Gomes Filho
