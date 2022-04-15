#!/bin/env python

import warnings
import os
import subprocess

################### Helper Functions ###########################################

def getPaired(sample, read, sample_dir):

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

def input4filter(wildcards, metadata):

    layout = metadata[wildcards.samples]["layout"]

    if layout.lower() == "single":
        return "results/{method}/{samples}/{samples}_{build}_se.bed"
    elif layout.lower() == "paired":    
        return "results/{method}/{samples}/{samples}_{build}_pe.bed"

def input4inpFasta(wildcards, metadata, sampleList):

    for sample in sampleList:

        if metadata[sample]["simulation_input"] == wildcards.samples:

            layout = metadata[sample]["simulation_input_layout"]

    if layout.lower() == "single":
        return "results/input/{samples}/{samples}_{build}_se.bed"
    elif layout.lower() == "paired":    
        return "results/input/{samples}/{samples}_{build}_pe.bed"

def input4PCA(sampleList, metadata, build, method):

    inputList = []
    for sample in sampleList:

        layout = metadata[sample]["layout"]

        if layout.lower() == "single":
            inputList.append(f"results/{method}/{sample}/{sample}_{build}_se_sortedbyCoordinates.bam")
        elif layout.lower() == "paired":   
            inputList.append(f"results/{method}/{sample}/{sample}_{build}_pe_sortedbyCoordinates.bam")

    return inputList

def input4nucTable(method):

    if method == "XR":
        return "results/XR/{samples}/{samples}_{build}_lengthMode.fa"
    elif method == "DS":    
        return "results/DS/{samples}/{samples}_{build}_sorted_10.fa"

def getMotif(sample, product):

    if product.lower() in ["oxaliplatin", "cisplatin", "bpdedg"]: 
        return "'.{4}(g|G){2}.{4}'"
    
    elif product.lower() in ["64", "64pp", "(6-4)pp", "6-4pp", "cpd"]: 
        return "'.{4}(c|t|C|T){2}.{4}'"

    elif product.lower() == "na":
        return "'.{10}'"

def getDinuc(sample, product):
    
    if product.lower() in ["oxaliplatin", "cisplatin"]: 
        return "'GG'"
    
    elif product.lower() in ["64", "64pp", "(6-4)pp", "6-4pp", "cpd"]: 
        return "'CC','CT','TC','TT'"


#def getInput(sample, inputExist, inputList, inputIdx, sampleList, build):
def getInput(sample, metadata, build):

    if "simulation_input" in metadata[sample]:

        input_name = metadata[sample]["simulation_input"]    
        return f"results/input/{input_name}/{input_name}_{build}.fasta"
    
    else:
        return f"resources/ref_genomes/{build}/genome_{build}.ron" 

def lineNum(file):
    
    linenum = 0
    if os.path.exists(file):
        with open(file) as f:
            for line in f:
                linenum += 1

    warnMessage = (f"\n{file} file is either empty or does not exists!\n" + 
        "It is expected if this is a dry-run. The file will be produced " + 
        "after the execution.")

    if linenum == 0:
        warnings.warn(warnMessage)

    return linenum

def mappedReads(*files):

    lineNumber = 0
    for file in files:
        lineNumber += lineNum(str(file))

    return lineNumber

def allInput(method, build, sampleList, metadata):

    inputList = []
    for sample in sampleList:
        sample_dir = f"results/{method}/{sample}/" 

        layout = metadata[sample]["layout"]
        simulation = metadata[sample]["simulation_enabled"]

        if layout.lower() == "single":
            inputList.append(f"{sample_dir}{sample}.html")
        elif layout.lower() == "paired":
            inputList.append(f"{sample_dir}{sample}_1.html")
            inputList.append(f"{sample_dir}{sample}_2.html")

        inputList.append(f"{sample_dir}{sample}_{build}_sorted_nucleotideTable.pdf")
        inputList.append(f"{sample_dir}{sample}_{build}_sorted_dinucleotideTable.pdf")
        inputList.append(f"{sample_dir}{sample}_{build}_length_distribution.pdf")
        inputList.append(f"{sample_dir}{sample}_{build}_{method}_sorted_plus.bw")
        inputList.append(f"{sample_dir}{sample}_{build}_{method}_sorted_minus.bw")

        if simulation:
            inputList.append(f"{sample_dir}{sample}_{build}_{method}_sim_nucleotideTable.pdf")
            inputList.append(f"{sample_dir}{sample}_{build}_{method}_sim_dinucleotideTable.pdf")

        if method == "DS":
            inputList.append(f"{sample_dir}{sample}_{build}_sorted_ds_dipyrimidines_plus.bed") 
            inputList.append(f"{sample_dir}{sample}_{build}_sorted_ds_dipyrimidines_minus.bed") 

        elif method == "XR":
            inputList.append(f"{sample_dir}{sample}_{build}_sorted_plus.bed") 
            inputList.append(f"{sample_dir}{sample}_{build}_sorted_minus.bed") 
    
    if method == "DS":
        inputList.append("results/scatterplot_PearsonCorr_bigwigScores_DS.pdf")
        inputList.append("results/PearsonCorr_bigwigScores_DS.tab")
        inputList.append("results/heatmap_SpearmanCorr_readCounts_DS.pdf")
        inputList.append("results/SpearmanCorr_readCounts_DS.tab")
        inputList.append("results/PCA_readCounts_DS.pdf")
    elif method == "XR":
        inputList.append("results/scatterplot_PearsonCorr_bigwigScores_XR.pdf")
        inputList.append("results/PearsonCorr_bigwigScores_XR.tab")
        inputList.append("results/heatmap_SpearmanCorr_readCounts_XR.pdf")
        inputList.append("results/SpearmanCorr_readCounts_XR.tab")
        inputList.append("results/PCA_readCounts_XR.pdf")
        
    #print(inputList)
    return inputList

################################################################################
