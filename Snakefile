#!/bin/env python

#### configuration file ####
configfile: "config/config_initial.yaml" 
configfile: "config/config.yaml"

# singularity image to use
containerized: "docker://azgarian/snakemake:1.2"

include: "workflow/rules/common.smk"

wildcard_constraints:
    build=config["genome"]["build"],
    method="|".join(set(config["meta"][sample]["method"] for sample in config["meta"])),

ruleorder: rename_raw > sra_pe
ruleorder: rename_raw_input > sra_pe_input
ruleorder: rename_raw > sra_se
ruleorder: rename_raw_input > sra_se_input

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