#!/bin/env python

import snakemake
version = snakemake.__version__

#### configuration file ####
configfile: "config/config_initial.yaml" 
configfile: "config/config.yaml"

is_profile_used = snakemake.snakemake.options.profile is not None
if is_profile_used:
    # If Snakemake version is 8 or above
    if version >= '8.0':
        # Configure for Snakemake >= 8.0
        config['executor'] = "cluster-generic"
        config['cluster-generic-submit-cmd'] = "sbatch -A {resources.account} -p {resources.partition} -J {rule}.job --qos {resources.partition} --cpus-per-task={threads} -e logs/cluster/{rule}_%A.err --output=/dev/null"
        config['software-deployment-method'] = "conda"
    else:
        # Configure for Snakemake < 8.0
        config['cluster'] = "sbatch -A {resources.account} -p {resources.partition} -J {rule}.job --qos {resources.partition} --cpus-per-task={threads} -e logs/cluster/{rule}_%A.err --output=/dev/null"
        config['use-conda'] = True

# singularity image to use
containerized: "docker://azgarian/snakemake:1.2"

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