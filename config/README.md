# Parameters of config_initial.yaml file

You can find the initial config example and the description of each parameter 
below:

```
DS:
  adaptor_se: "-g GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT"
  adaptor_pe: "-g GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT -G GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT"
  cutadapt_se: "--discard-trimmed"
  cutadapt_pe: "--discard-trimmed"
  samtools_se: "-q 20"
  samtools_pe: "-q 20 -bf 0x2"

XR:
  adaptor_se: "-a TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG"
  adaptor_pe: "-a TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG -A TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG"
  cutadapt_se: ""
  cutadapt_pe: ""
  samtools_se: "-q 20"
  samtools_pe: "-q 20 -bf 0x2"
  
#######

filter: "'^'"
```

- `DS`: contains the parameter settings unique to Damage-seq analysis.

- `XR`: contains the parameter settings unique to XR-seq analysis.

- `filter`: A regex that will be used to filter the aligned reads 
    (in BED format). If no filtering will be done, `"'^'"` regex can be used 
    to grep all the lines. 

All the parameters below are used for both `DS` and `XR`:

    - `adaptor_se`: single-end adaptor parameters of cutadapt (-a, -g...).
    - `adaptor_pe`: paired-end adaptor parameters of cutadapt (-a, -g, -A, -G...).
    -  cutadapt_se: any extra parameters of cutadapt for single-end samples.
    -  cutadapt_pe: any extra parameters of cutadapt for paired-end samples.
    -  samtools_se: samtools flags for quality filtering (single-end).
    -  samtools_pe: samtools flags for quality filtering (paired-end).