#!/bin/env python

#### configuration file ####
configfile: "config/config_XR_initial.yaml" 
configfile: "config/config_XR.yaml"

include: "workflow/rules/common.smk"

rule all:
    input:
        lambda w: allInput(config["method"], config["build"], config["sample"], 
            config["meta"]),

if config["srr"]["enabled"]:
    include: "workflow/rules/sra.smk"
else:    
    include: "workflow/rules/rename_raw.smk"

if config["input"]["srr"]["enabled"]:
    include: "workflow/rules/sra_input.smk"
else:
    include: "workflow/rules/rename_raw_input.smk"
    
include: "workflow/rules/fastqc.smk"

if config["genome_download"]:
    include: "workflow/rules/genome_download.smk"

if config["bowtie2_build"]:
    include: "workflow/rules/genome_build.smk"
    
include: "workflow/rules/genome_indexing.smk"
include: "workflow/rules/genome_idx2ron.smk"
include: "workflow/rules/adaptor_handling.smk"
include: "workflow/rules/align.smk"
include: "workflow/rules/bam2bed.smk"
include: "workflow/rules/sort_filter.smk"
include: "workflow/rules/length_distribution.smk"
include: "workflow/rules/plot_length.smk"
include: "workflow/rules/length_mode.smk"
include: "workflow/rules/sep_strands.smk"
include: "workflow/rules/bed2fasta.smk"
include: "workflow/rules/nucleotide_table.smk"
include: "workflow/rules/plot_nuc.smk"
include: "workflow/rules/genomecov.smk"
include: "workflow/rules/bedGraphToBigWig.smk"
include: "workflow/rules/simulation.smk"
include: "workflow/rules/nucleotide_table_sim.smk"
include: "workflow/rules/plot_nuc_sim.smk"
#include: "workflow/rules/bgzip_variation.smk"
#include: "workflow/rules/igv_reports.smk"
include: "workflow/rules/bam_correlation.smk"
include: "workflow/rules/bam_corr_graphs.smk"


