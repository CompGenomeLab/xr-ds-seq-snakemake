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
        #input_list.append(f"results/XR/{sample}/{sample}_cutadapt_se_hg38.bam")
        input_list.append(f"results/XR/{sample}/{sample}_T2C_cutadapt_se_hg38.bam")
        input_list.append(f"results/XR/{sample}/{sample}_masked_cutadapt_se_hg38.bam")
        # input_list.append(f"results/XR/{sample}/{sample}_hg38_se_sortedbyCoordinates.bam")
        # input_list.append(f"results/XR/{sample}/{sample}_cutadapt_se_very_sensitive_hg38.bam")
        # input_list.append(f"results/XR/{sample}/{sample}_cutadapt_se_bwa_hg38.bam")
        # input_list.append(f"results/XR/{sample}/{sample}_cutadapt_se_bwa_sensitive1_hg38.bam")
        # input_list.append(f"results/XR/{sample}/{sample}_cutadapt_se_bwa_sensitive2_hg38.bam")
        # input_list.append(f"results/XR/{sample}/{sample}_cutadapt_se_bwa_sensitive3_hg38.bam")
        # input_list.append(f"results/XR/{sample}/{sample}_cutadapt_se_bwa_sensitive4_hg38.bam")
        # input_list.append(f"results/XR/{sample}/{sample}_cutadapt_se_bwa_sensitive5_hg38.bam")
        # input_list.append(f"results/processed_files/{sample}_fastqc.zip")
        # input_list.append(f"results/processed_files/{sample}_cutadapt_fastqc.zip")
        # input_list.append(f"results/processed_files/{sample}_hg38_XR_dinucleotideTable.pdf")
    return input_list

samples = [
    "NHF1UV2hCPD_CTTATCGG_S12_R1_001",
    "NHF1UV4hCPD_TCCGCTAA_S13_R1_001",
    "NHF1UV8hCPD_GATCTATC_S14_R1_001",
    "NDDB2CPD2h_AGCGCTAG_S10_R1_001",
    "NDDB2CPD4h_GATATCGA_S12_R1_001",
    "NHF1CPD2h_GATCTATC_S14_R1_001",
    "NHF1CPD4h_AGCTCGCT_S16_R1_001",
    "NDDB2CPD30m_GTGTCGGA_S11_R1_001",
    "NHF1CPD10m_CTTATCGG_S13_R1_001",
    "NHF1CPD30m_TCCGCTAA_S15_R1_001",
    "NHF1CPD8h_ACACTAAG_S17_R1_001"
    ]

rule all:
    input:
        #lambda w: allInput(config["genome"]["build"], config["meta"]),
        lambda w: input_rule(samples),

# Downloading reference genome and generating related files  
include: "workflow/rules/prepare_genome.smk"

# Download and rename (if necessary) sample files
include: "workflow/rules/download_sample.smk"

# Quality control
include: "workflow/rules/fastqc.smk"

# Adaptor handling, mapping, quality trimming, and converting to bed
include: "workflow/rules/fastq2bed.smk"
#include: "workflow/rules/alignment.smk"

# Sorting, filtering, calculating length dist., filtering damage-seq samples by
# motif, producing bigwig files
include: "workflow/rules/process_bed.smk"

# Simulating the sample reads
include: "workflow/rules/simulation.smk"

# Plotting length distribution, nucleotide enrichment, and bam correlations
include: "workflow/rules/plot.smk"