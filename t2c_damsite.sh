awk 'NR % 4 == 1 {header = $0} 
  NR % 4 == 2 {
      seq = $0;
      changed = 0;
      if (substr(seq, length(seq)-6, 1) == "T") {
          seq = substr(seq, 1, length(seq)-7) "C" substr(seq, length(seq)-5);
          changed = 1;
      }
      if (substr(seq, length(seq)-7, 1) == "T") {
          seq = substr(seq, 1, length(seq)-8) "C" substr(seq, length(seq)-6);
          changed = 1;
      }
      if (changed) {
          print header;
          print seq;
          getline; print $0;
          getline; print $0;
      } else {
          getline; getline; 
      }
  }'  NHF1_CPD_1h_XR_rep1_mapped_as0.fastq > NHF1_CPD_1h_XR_rep1_mapped_as0_T2C_damsite.fastq

awk '{
    if (NR % 4 == 1) header = $0;         # Capture the header line
    else if (NR % 4 == 2) seq = $0;      # Capture the sequence line
    else if (NR % 4 == 0) qual = $0;     # Capture the quality line
    else if (NR % 4 == 3) plus = $0;     # Capture the "+" line
    if (NR % 4 == 0 && length(seq) == length(qual)) {
        print header; print seq; print plus; print qual; # Output valid reads
    }
}' NHF1_CPD_1h_XR_rep1_mapped_as0_T2C_damsite.fastq > NHF1_CPD_1h_XR_rep1_mapped_as0_T2C_damsite2.fastq

mv NHF1_CPD_1h_XR_rep1_mapped_as0_T2C_damsite2.fastq NHF1_CPD_1h_XR_rep1_mapped_as0_T2C_damsite.fastq

gzip NHF1_CPD_1h_XR_rep1_mapped_as0_T2C_damsite.fastq