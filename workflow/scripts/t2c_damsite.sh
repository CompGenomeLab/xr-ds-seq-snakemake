#!/bin/bash

# Usage: ./process_fastq.sh <input.fastq> <output.fastq.gz>

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.fastq> <output.fastq.gz>"
    exit 1
fi

input_fastq="$1"
output_fastq="$2"
temp_fastq="${output_fastq%.gz}"  # Remove .gz if present

awk 'NR % 4 == 1 {header = $0} 
  NR % 4 == 2 {
      seq = $0;
      if (substr(seq, length(seq)-6, 1) == "T") {
          seq = substr(seq, 1, length(seq)-7) "C" substr(seq, length(seq)-5);
      }
      if (substr(seq, length(seq)-7, 1) == "T") {
          seq = substr(seq, 1, length(seq)-8) "C" substr(seq, length(seq)-6);
      }
      print header;
      print seq;
      getline; print $0;
      getline; print $0;
  }'  "$input_fastq" > "$temp_fastq"

awk '{
    if (NR % 4 == 1) header = $0;         # Capture the header line
    else if (NR % 4 == 2) seq = $0;      # Capture the sequence line
    else if (NR % 4 == 0) qual = $0;     # Capture the quality line
    else if (NR % 4 == 3) plus = $0;     # Capture the "+" line
    if (NR % 4 == 0 && length(seq) == length(qual)) {
        print header; print seq; print plus; print qual; # Output valid reads
    }
}' "$temp_fastq" > "${temp_fastq}.tmp"

mv "${temp_fastq}.tmp" "$temp_fastq"

# Step 3: Gzip if output ends with .gz
if [[ "$output_fastq" == *.gz ]]; then
    gzip -c "$temp_fastq" > "$output_fastq"
    rm "$temp_fastq"
else
    mv "$temp_fastq" "$output_fastq"
fi

echo "Processed FASTQ saved to $output_fastq"