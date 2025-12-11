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
        # input_list.append(f"results/XR/{sample}/{sample}_masked_cutadapt_se_hg38.bam")
        # input_list.append(f"results/XR/{sample}/{sample}_hg38_se_sortedbyCoordinates.bam")
        input_list.append(f"results/XR/{sample}/{sample}_hg38.bam")
        # input_list.append(f"results/processed_files/{sample}_fastqc.zip")
        # input_list.append(f"results/processed_files/{sample}_cutadapt_fastqc.zip")
        # input_list.append(f"results/XR/{sample}/{sample}_cutadapt_se_hg38_unmapped_filtered_overrep_fastqc.zip")
        # input_list.append(f"results/processed_files/{sample}_hg38_XR_filtered_dinucleotideTable.pdf")
    return input_list


samples = [
"CSBUVB2hHirtTaq_GATATCGA_S1",
"CSBUVC2hHirtTaq_AGCTCGCT_S2",
"DDB2UVB2hHirtTaq_GTGTCGGA_S3",  
"DDB2UVC2hHirtTaq_TCCGCTAA_S4",  
"NHF1UVB2hHirtTaq_ACACTAAG_S5",
"NHF1UVC2hHirtTaq_CTTATCGG_S6",
"XPCUVB2hHirtTaq_AGCGCTAG_S19",
"XPCUVC2hHirtTaq_GATCTATC_S20",
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
include: "workflow/rules/fastq2bed_alt.smk"
#include: "workflow/rules/alignment.smk"

# Sorting, filtering, calculating length dist., filtering damage-seq samples by
# motif, producing bigwig files
include: "workflow/rules/process_bed.smk"

# Simulating the sample reads
include: "workflow/rules/simulation.smk"

# Plotting length distribution, nucleotide enrichment, and bam correlations
include: "workflow/rules/plot.smk"