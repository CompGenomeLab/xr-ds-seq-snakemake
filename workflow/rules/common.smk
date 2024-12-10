#!/bin/env python

import warnings
import os
import subprocess

################### Helper Functions ###########################################

def getMethodParams(wildcards, metadata, parameter, XR, DS):

    """
    Retrieves the given parameters from config_initial.yaml. 

    Used rules: cutadapt_se, cutadapt_pe, bam2bed_se, bam2bed_pe
    """

    method = metadata[wildcards.samples]["method"].upper()
    layout = metadata[wildcards.samples]["layout"].lower()
    lo = "se" if layout == "single" else "pe"

    return XR[f"{parameter.lower()}_{lo}"] if method == "XR" else DS[f"{parameter.lower()}_{lo}"]

def getPossibleInputNames(metadata):

    """
    Gets all the possible SRR ids for input files if there is one. 

    Used rules: sra_se_input, sra_pe_input
    """  

    input_list = []
    for sample in metadata.keys():
        if "input" in metadata[sample]["simulation"]:
            input_list.append(metadata[sample]["simulation"]["input"]["name"])

    return input_list

def getSRR(wildcards, metadata):

    """
    Gets the SRR id if there is one. 

    Used rules: sra_se, sra_pe
    """
    
    if "srr_id" in metadata[wildcards.samples]:
        return metadata[wildcards.samples]["srr_id"]
    else:
        return " "

def input4rule(wildcards, metadata, rule, filtered=False, build=config["genome"]["build"]):

    """
    Retreives the input file/s for given rule.

    Used rules: sort_filter, bed2fasta_input, bam_correlation, nucleotide_table,
                simulation_ds, simulation_xr
    """

    if rule == "sort_filter":

        layout = metadata[wildcards.samples]["layout"].lower()
        lo = "se" if layout == "single" else "pe"
        return "results/{method}/{samples}/{samples}_{build}" + f"_{lo}.bed"

    elif rule == "bed2fasta_input":

        for sample, value in metadata.items():
            if "input" in value["simulation"]:
                if value["simulation"]["input"]["name"] == wildcards.samples:
                    layout = value["simulation"]["input"]["layout"].lower()
                    break

        lo = "se" if layout == "single" else "pe"
        return "results/input/{samples}/{samples}_{build}" + f"_{lo}.bed"

    elif rule == "bam_correlation":

        inputList = []
        for sample in metadata.keys():
            lo = "se" if metadata[sample]["layout"].lower() == "single" else "pe"
            method = metadata[sample]["method"].upper()
            inputList.append(f"results/{method}/{sample}/{sample}_cutadapt_{lo}_{build}_sorted.bam")

        return inputList

    elif rule == "bam_correlation_idx":

        inputList = []
        for sample in metadata.keys():
            lo = "se" if metadata[sample]["layout"].lower() == "single" else "pe"
            method = metadata[sample]["method"].upper()
            inputList.append(f"results/{method}/{sample}/{sample}_cutadapt_{lo}_{build}_sorted.bam.bai")

        return inputList

    elif rule == "nucleotide_table":

        method = metadata[wildcards.samples]["method"].upper()
        if method == "XR":
            ext = "_lengthMode.fa"
        elif method == "DS":
            ext = "_sorted_filt_10.fa" if filtered else "_sorted_10.fa"

        return "results/{method}/{samples}/{samples}_{build}" + ext

    elif rule == "simulation_xr" or rule == "simulation_ds":

        if "input" in metadata[wildcards.samples]["simulation"]:
            input_name = metadata[wildcards.samples]["simulation"]["input"]["name"]    
            return f"results/input/{input_name}/{input_name}_{build}.fasta"      
        else:
            return f"resources/ref_genomes/{build}/genome_{build}.ron" 

def getMotif(sample, product):

    """
    Filters reads with a motif based on the damage type of the damage-seq 
    sample.

    Used rules: filtbyMotifs
    """

    if product.lower().strip() in ["oxaliplatin", "cisplatin", "bpdedg"]: 
        return "'.{4}(g|G){2}.{4}'"
    
    elif product.lower().strip() in ["64", "64pp", "(6-4)pp", "6-4pp", "64-pp", "cpd"]: 
        return "'.{4}(c|t|C|T){2}.{4}'"

    elif product.lower().strip() == "na":
        return "'.{10}'"

def getDinuc(sample, product):
    
    """
    Selects the nucleotides to plot based on the damage type of the sample.

    Used rules: plot_nuc, plot_nuc_sim
    """

    if product.lower() in ["oxaliplatin", "cisplatin"]: 
        return "'GG'"
    
    elif product.lower() in ["64", "64pp", "(6-4)pp", "6-4pp", "cpd"]: 
        return "'CC','CT','TC','TT'"

def lineNum(file):

    """
    Calculates the line number of a given file.
    """
    
    linenum = 0
    if os.path.exists(file):
        with open(file) as f:
            for line in f:
                linenum += 1

    warnMessage = (f"\n{file} file is empty!\n")

    if linenum == 0 and os.path.exists(file):
        warnings.warn(warnMessage)

    return linenum

def mappedReads(*files):

    """
    Sums the line number of all given files. Used for calculated mapped reads
    of a sample.

    Used rules: genomecov_ds, genomecov_xr, genomecov_sim
    """

    lineNumber = 0
    for file in files:
        lineNumber += lineNum(str(file))

    return lineNumber

def allInput(build, metadata):

    """
    Returns the names of the output files to be produced when the pipeline ends.

    Used rules: all (this is the name of the rule)
    """

    inputList = []
    for sample in metadata.keys():

        method = metadata[sample]["method"].upper()
        sdir = "results/processed_files" 
        sprefix = f"{sample}_{build}_{method}"
        layout = metadata[sample]["layout"]
        simulation = metadata[sample]["simulation"]["enabled"]


        if layout.lower() == "single":
            inputList.append(f"{sdir}/{sample}_fastqc.html")
        elif layout.lower() == "paired":
            inputList.append(f"{sdir}/{sample}_1_fastqc.html")
            inputList.append(f"{sdir}/{sample}_2_fastqc.html")

        inputList.append(f"{sdir}/{sprefix}_nucleotideTable.pdf")
        inputList.append(f"{sdir}/{sprefix}_dinucleotideTable.pdf")
        inputList.append(f"{sdir}/{sprefix}_filt_nucleotideTable.pdf")
        inputList.append(f"{sdir}/{sprefix}_filt_dinucleotideTable.pdf")
        inputList.append(f"{sdir}/{sprefix}_length_distribution.pdf")
        inputList.append(f"{sdir}/{sprefix}_plus.bed") 
        inputList.append(f"{sdir}/{sprefix}_minus.bed") 
        inputList.append(f"{sdir}/{sprefix}_plus.bw")
        inputList.append(f"{sdir}/{sprefix}_minus.bw")

        if simulation:
            inputList.append(f"{sdir}/{sprefix}_sim_nucleotideTable.pdf")
            inputList.append(f"{sdir}/{sprefix}_sim_dinucleotideTable.pdf")
            inputList.append(f"{sdir}/{sprefix}_sim_plus.bed")
            inputList.append(f"{sdir}/{sprefix}_sim_minus.bed")
            inputList.append(f"{sdir}/{sprefix}_sim_plus.bw")
            inputList.append(f"{sdir}/{sprefix}_sim_minus.bw")
    
    if len(metadata.keys()) > 1:

        inputList.append(f"{sdir}/scatterplot_PearsonCorr_bigwigScores.pdf")
        inputList.append(f"{sdir}/PearsonCorr_bigwigScores.tab")
        inputList.append(f"{sdir}/heatmap_SpearmanCorr_readCounts.pdf")
        inputList.append(f"{sdir}/SpearmanCorr_readCounts.tab")
        inputList.append(f"{sdir}/PCA_readCounts.pdf")
            
    #print(inputList)
    return inputList

################################################################################
