#!/bin/env python

import warnings
import os
import subprocess

################### Helper Functions ###########################################

def getPaired(sample, read, sample_dir):

    """
    Finds the suffix of a given sample name.

    Used rules: rename_raw, rename_raw_input 
    """

    pairedR1 = f"{sample_dir}{sample}_R1.fastq.gz"
    paired1 = f"{sample_dir}{sample}_1.fastq.gz"
    
    if os.path.isfile(pairedR1) and read == "forward":
        return f"{sample_dir}{sample}_R1.fastq.gz"

    elif os.path.isfile(pairedR1) and read == "reverse":
        return f"{sample_dir}{sample}_R2.fastq.gz"

    elif os.path.isfile(paired1) and read == "forward":
        return f"{sample_dir}{sample}_1.fastq.gz"

    elif os.path.isfile(paired1) and read == "reverse":
        return f"{sample_dir}{sample}_2.fastq.gz"

def getMethodParams(wildcards, metadata, parameter, XR, DS):

    """
    Retrieves the given parameters from config_initial.yaml. 

    Used rules: cutadapt_se, cutadapt_pe, bam2bed_se, bam2bed_pe
    """

    method = metadata[wildcards.samples]["method"].upper()
    layout = metadata[wildcards.samples]["layout"].lower()

    if parameter.lower() == "adaptor" and layout== "single":
        if method == "XR":
            return XR["adaptor_se"]
        elif method == "DS":
            return DS["adaptor_se"]
    elif parameter.lower() == "adaptor" and layout== "paired":
        if method == "XR":
            return XR["adaptor_pe"]
        elif method == "DS":
            return DS["adaptor_pe"]
    elif parameter.lower() == "cutadapt" and layout== "single":
        if method == "XR":
            return XR["cutadapt_se"]
        elif method == "DS":
            return DS["cutadapt_se"]
    elif parameter.lower() == "cutadapt" and layout== "paired":
        if method == "XR":
            return XR["cutadapt_pe"]
        elif method == "DS":
            return DS["cutadapt_pe"]
    elif parameter.lower() == "samtools" and layout== "single":
        if method == "XR":
            return XR["samtools_se"]
        elif method == "DS":
            return DS["samtools_se"]
    elif parameter.lower() == "samtools" and layout== "paired":
        if method == "XR":
            return XR["samtools_pe"]
        elif method == "DS":
            return DS["samtools_pe"]

def getSRR(wildcards, metadata):

    """
    Gets the SRR id if there is one. 

    Used rules: sra_se, sra_pe
    """

    if "srr_id" in metadata[wildcards.samples]:
        return metadata[wildcards.samples]["srr_id"]
    else:
        return " "

def input4filter(wildcards, metadata):

    """
    Retreives the input file of sort_filter rule.

    Used rules: sort_filter
    """

    layout = metadata[wildcards.samples]["layout"].lower()

    if layout == "single":
        return "results/{method}/{samples}/{samples}_{build}_se.bed"
    elif layout == "paired":    
        return "results/{method}/{samples}/{samples}_{build}_pe.bed"

def input4inpFasta(wildcards, metadata):

    """
    Retreives the input file of bed2fasta_input rule based on the simulation 
    info in the config file.

    Used rules: bed2fasta_input
    """

    for sample in metadata.keys():

        if metadata[sample]["simulation"]["input"]["name"] == wildcards.samples:
            layout = metadata[sample]["simulation"]["input"]["layout"]

    if layout.lower() == "single":
        return "results/input/{samples}/{samples}_{build}_se.bed"
    elif layout.lower() == "paired":    
        return "results/input/{samples}/{samples}_{build}_pe.bed"

def input4PCA(metadata, build):

    """
    Retreives the input file of bam_correlation rule.

    Used rules: bam_correlation
    """

    inputList = []
    for sample in metadata.keys():

        layout = metadata[sample]["layout"].lower()
        method = metadata[sample]["method"].upper()

        if layout == "single":
            inputList.append(f"results/{method}/{sample}/{sample}_{build}_se_sortedbyCoordinates.bam")
        elif layout == "paired":   
            inputList.append(f"results/{method}/{sample}/{sample}_{build}_pe_sortedbyCoordinates.bam")

    return inputList

def input4nucTable(wildcards, metadata):

    """
    Retreives the input file of nucleotide_table rule.

    Used rules: nucleotide_table
    """

    method = metadata[wildcards.samples]["method"].upper()

    if method == "XR":
        return "results/XR/{samples}/{samples}_{build}_lengthMode.fa"
    elif method == "DS":    
        return "results/DS/{samples}/{samples}_{build}_sorted_10.fa"

def getMotif(sample, product):

    """
    Filters reads with a motif based on the damage type of the damage-seq 
    sample.

    Used rules: filtbyMotifs
    """

    if product.lower() in ["oxaliplatin", "cisplatin", "bpdedg"]: 
        return "'.{4}(g|G){2}.{4}'"
    
    elif product.lower() in ["64", "64pp", "(6-4)pp", "6-4pp", "cpd"]: 
        return "'.{4}(c|t|C|T){2}.{4}'"

    elif product.lower() == "na":
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

def getInput(sample, metadata, build):

    """
    Retrieves the input file of the sample for simulation.

    Used rules: simulation_ds, simulation_xr
    """

    if "input" in metadata[sample]["simulation"]:

        input_name = metadata[sample]["simulation"]["input"]["name"]    
        return f"results/input/{input_name}/{input_name}_{build}.fasta"
    
    else:
        return f"resources/ref_genomes/{build}/genome_{build}.ron" 

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
