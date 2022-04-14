#!/bin/env python

import warnings
import os
import subprocess

################### Helper Functions ###########################################

def isSingle(sample, srrEnabled, srrID, sample_dir):

    if srrEnabled:

        if srrID == "NA":

            single = f"{sample_dir}{sample}.fastq.gz"
            pairedR1 = f"{sample_dir}{sample}_R1.fastq.gz"
            paired1 = f"{sample_dir}{sample}_1.fastq.gz"
            
            if os.path.isfile(pairedR1) or os.path.isfile(paired1):
                return False
            elif os.path.isfile(single):
                return True
            else:
                raise(ValueError(f"{paired1}, {pairedR1}, or {single} not found..."))

        if ":" in srrID:
            srrID = srrID.split(":")[0]

        shellCommand = f'fastq-dump -X 1 -Z --split-spot {srrID} | wc -l'
        #print(shellCommand)
        p=subprocess.getoutput(shellCommand)
        #print(p)
        lineNum = int(p.split("\n")[2])
        #print(lineNum)

        if lineNum == 4:
            return True
        else:
            return False

    else:

        single = f"{sample_dir}{sample}.fastq.gz"
        pairedR1 = f"{sample_dir}{sample}_R1.fastq.gz"
        paired1 = f"{sample_dir}{sample}_1.fastq.gz"
        
        if os.path.isfile(pairedR1) or os.path.isfile(paired1):
            return False
        elif os.path.isfile(single):
            return True
        else:
            raise(ValueError(f"{paired1}, {pairedR1}, or {single} not found..."))

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

    srrEnabled = metadata[wildcards.samples]["srr_enabled"]
    srrID = metadata[wildcards.samples]["srr_id"]

    if isSingle(wildcards.samples, srrEnabled, srrID, "resources/samples/"):
        return "results/{method}/{samples}/{samples}_{build}_se.bed"
    else:    
        return "results/{method}/{samples}/{samples}_{build}_pe.bed"

def input4inpFasta(wildcards, metadata):

    srrEnabled = metadata[wildcards.samples]["srr_enabled"]
    srrID = metadata[wildcards.samples]["srr_id"]

    if isSingle(wildcards.samples, srrEnabled, srrID, "resources/input/"):
        return "results/input/{samples}/{samples}_{build}_se.bed"
    else:    
        return "results/input/{samples}/{samples}_{build}_pe.bed"

def input4PCA(sampleList, metadata, build, method):

    inputList = []
    for sample in sampleList:

        srrEnabled = metadata[sample]["srr_enabled"]
        srrID = metadata[sample]["srr_id"]

        if isSingle(sample, srrEnabled, srrID, "resources/samples/"):
            inputList.append(f"results/{method}/{sample}/{sample}_{build}_se_sortedbyCoordinates.bam")
        else:    
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

def getInput(sample, inputExist, inputList, inputIdx, sampleList, build):

    if inputExist:
        inpDict={}
        for inp_idx in range(len(inputIdx)):
            idx_split = inputIdx[inp_idx].strip().split(",")
            indexList=[]
            for sample_idx in idx_split:
                sample_idx = sample_idx.strip() 
                if "-" in sample_idx:
                    for range_idx in range(int(sample_idx.split("-")[0]), int(sample_idx.split("-")[1])+1):
                        indexList.append(int(range_idx)) 
                else:
                    indexList.append(int(sample_idx)) 
                
            for sample_idx in indexList:
                if inputList[inp_idx] not in inpDict:
                    inpDict[inputList[inp_idx]] = [sampleList[sample_idx]]
                else:
                    inpDict[inputList[inp_idx]].append(sampleList[sample_idx])
        for k,v in inpDict.items():
        
            if sample in v:
            
                return f"results/input/{k}/{k}_{build}.fasta"
                
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

        srrEnabled = metadata[sample]["srr_enabled"]
        srrID = metadata[sample]["srr_id"]

        if isSingle(sample, srrEnabled, srrID, "resources/samples/"):
            inputList.append(f"{sample_dir}{sample}.html")
        else:
            inputList.append(f"{sample_dir}{sample}_1.html")
            inputList.append(f"{sample_dir}{sample}_2.html")

        inputList.append(f"{sample_dir}{sample}_{build}_sorted_nucleotideTable.pdf")
        inputList.append(f"{sample_dir}{sample}_{build}_sorted_dinucleotideTable.pdf")
        inputList.append(f"{sample_dir}{sample}_{build}_{method}_sim_nucleotideTable.pdf")
        inputList.append(f"{sample_dir}{sample}_{build}_{method}_sim_dinucleotideTable.pdf")
        inputList.append(f"{sample_dir}{sample}_{build}_length_distribution.pdf")
        inputList.append(f"{sample_dir}{sample}_{build}_{method}_sorted_plus.bw")
        inputList.append(f"{sample_dir}{sample}_{build}_{method}_sorted_minus.bw")
        #inputList.append(f"{sample_dir}{sample}_{build}_igv_report_chrX.html")

        #for i in range(1,23):
        #    inputList.append(f"{sample_dir}{sample}_{build}_igv_report_chr{str(i)}.html")

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
