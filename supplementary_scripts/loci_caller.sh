#!/bin/bash

# Load TTS sites into an associative array
declare -A TTS

# Read each TTS entry into the array
while read -r chr source type start end score strand phase _; do
    if [[ "$type" == "TTS" ]]; then
        key="$chr:$strand"
        TTS["$key"]+="$start,$end;"
    fi
done < tts.gff

# Process TSS sites and find matching TTS sites within 400 bp on the same strand
while read -r chr source type start end score strand phase attributes; do
    if [[ "$type" == "TSS" ]]; then
        key="$chr:$strand"
        if [[ -n "${TTS[$key]}" ]]; then
            max_distance=150
            furthest_tts_end=0
            furthest_distance=0
            IFS=';' read -ra tts_sites <<< "${TTS[$key]}"
            for tts_site in "${tts_sites[@]}"; do
                IFS=',' read -r tts_start tts_end <<< "$tts_site"
                
                if [[ "$strand" == "+" ]]; then
                    # For the "+" strand, TTS must be downstream of TSS
                    distance=$((tts_start - start))
                    if ((distance > 0 && distance <= max_distance)); then
                        if ((distance > furthest_distance)); then
                            furthest_distance=$distance
                            furthest_tts_end=$tts_end
                        fi
                    fi
                elif [[ "$strand" == "-" ]]; then
                    # For the "-" strand, TTS must be upstream of TSS
                    distance=$((start - tts_end))
                    if ((distance > 0 && distance <= max_distance)); then
                        if ((distance > furthest_distance)); then
                            furthest_distance=$distance
                            furthest_tts_end=$tts_start
                        fi
                    fi
                fi
            done

            # Output result if a valid TTS was found within 400 bp
            if ((furthest_tts_end > 0)); then
                if [[ "$strand" == "+" ]]; then
                    # For the "+" strand, output the transcript range from TSS to TTS end
                    echo -e "$chr\t$source\ttranscript\t$start\t$furthest_tts_end\t.\t$strand\t.\tTSS=$attributes;TTS_end=$furthest_tts_end"
                elif [[ "$strand" == "-" ]]; then
                    # For the "-" strand, output the transcript range from TTS start to TSS
                    echo -e "$chr\t$source\ttranscript\t$furthest_tts_end\t$start\t.\t$strand\t.\tTSS=$attributes;TTS_start=$furthest_tts_end"
                fi
            fi
        fi
    fi
done < tss.gff

#post processing
#bedtools merge -s -i t1.gff -c 6,7,8 -o distinct
#adding asRNA tag
#awk '{print $1 "\tTSSAR\tgene\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\tID=asRNA_" NR ";Name=asRNA_" NR ";locus_tag=asRNA_" NR}' input_file > output_file
#awk '{print $0 "\tID=asRNA_" NR ";Name=asRNA_" NR ";locus_tag=asRNA_" NR}' input_file > output_file
