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
# motif
include: "workflow/rules/process_bed.smk"

# Simulating the sample reads
include: "workflow/rules/simulation.smk"

# BigWig files for igv
include: "workflow/rules/bed2bigWig.smk"

# Plots
include: "workflow/rules/plot_length.smk"
include: "workflow/rules/plot_nuc.smk"
include: "workflow/rules/plot_bam_corr.smk"

# Copy bed files to final destination
include: "workflow/rules/copy_bed.smk"