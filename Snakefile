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

def input_rule(samplist):
    input_list = []
    for sample in samplist:
        input_list.append(f"results/XR/{sample}/{sample}_cutadapt_se_hg38.bam")
        input_list.append(f"results/XR/{sample}/{sample}_cutadapt_se_very_sensitive_hg38.bam")
        input_list.append(f"results/XR/{sample}/{sample}_cutadapt_se_bwa_hg38.bam")
        input_list.append(f"results/XR/{sample}/{sample}_cutadapt_se_bwa_sensitive1_hg38.bam")
        input_list.append(f"results/XR/{sample}/{sample}_cutadapt_se_bwa_sensitive2_hg38.bam")
        input_list.append(f"results/XR/{sample}/{sample}_cutadapt_se_bwa_sensitive3_hg38.bam")
        input_list.append(f"results/XR/{sample}/{sample}_cutadapt_se_bwa_sensitive4_hg38.bam")
        input_list.append(f"results/XR/{sample}/{sample}_cutadapt_se_bwa_sensitive5_hg38.bam")
    return input_list

samples = [
    "NHF1_CPD_1h_XR_rep1",
    "NHF1_64_1h_XR_rep1",
    "NHF1UV2h64Kappa_TCCGCTAA_S16_L001_R1_001",
    "NHF1UV2hCPDKappa_GATCTATC_S19_L001_R1_001",
    ]

rule all:
    input:
        lambda w: allInput(config["genome"]["build"], config["meta"]),
        #lambda w: input_rule(samples),

# Downloading reference genome and generating related files  
include: "workflow/rules/prepare_genome.smk"

# Download and rename (if necessary) sample files
include: "workflow/rules/download_sample.smk"

# Quality control
include: "workflow/rules/fastqc.smk"

# Adaptor handling, mapping, quality trimming, and converting to bed
include: "workflow/rules/fastq2bed.smk"
# include: "workflow/rules/alignment.smk"

# Sorting, filtering, calculating length dist., filtering damage-seq samples by
# motif, producing bigwig files
include: "workflow/rules/process_bed.smk"

# Simulating the sample reads
include: "workflow/rules/simulation.smk"

# Plotting length distribution, nucleotide enrichment, and bam correlations
include: "workflow/rules/plot.smk"