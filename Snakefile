#!/bin/env python

#### configuration file ####
configfile: "config/config_initial.yaml" 
configfile: "config/config.yaml"

include: "workflow/rules/common.smk"

wildcard_constraints:
    build=config["genome"]["build"],
    method="XR|DS"
    
rule all:
    input:
        lambda w: allInput(config["genome"]["build"], config["meta"]),

# Downloading reference genome and generating related files  
include: "workflow/rules/prepare_genome.smk"

# Download and rename (if necessary) sample files
include: "workflow/rules/download_sample.smk"

# Quality control
include: "workflow/rules/fastqc.smk"

# Adaptor handling, mapping, quality trimming, and converting to bed
include: "workflow/rules/fastq2bed.smk"

# Sorting, filtering, calculating length dist., filtering damage-seq samples by
# motif, producing bigwig files
include: "workflow/rules/process_bed.smk"

# Simulating the sample reads
include: "workflow/rules/simulation.smk"

# Plotting length distribution, nucleotide enrichment, and bam correlations
include: "workflow/rules/plot.smk"