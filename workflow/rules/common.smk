#!/bin/env python

import warnings
import os
import subprocess

################### Helper Functions ###########################################

def getSRR(sample, srrList, sampleList):

    try:
        idx = sampleList.index(sample)
    except:
       raise(ValueError("Designated wildcard cannot be found in sample list."))
        
    return srrList[idx]

def isSingle(sample, sampleList, srrEnabled, srrList):

    if srrEnabled:

        mySRR = getSRR(sample, srrList, sampleList)

        if ":" in mySRR:
            mySRR = mySRR.split(":")[0]

        shellCommand = 'fastq-dump -X 1 -Z --split-spot ' + mySRR + ' | wc -l'
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

        sample_dir = "resources/samples/"
        single = sample_dir + sample + ".fastq.gz"
        pairedR1 = sample_dir + sample + "_R1.fastq.gz"
        paired1 = sample_dir + sample + "_1.fastq.gz"
        
        if os.path.isfile(pairedR1) or os.path.isfile(paired1):
            return False
        elif os.path.isfile(single):
            return True
        else:
            raise(ValueError(pairedR1, single, "Sample not found..."))

def getPaired(sample, sampleList, read):

    sample_dir = "resources/samples/"
    pairedR1 = sample_dir + sample + "_R1.fastq.gz"
    #paired1 = sample_dir + sample + "_1.fastq.gz"
    
    if os.path.isfile(pairedR1) and read == "forward":
        return sample_dir + sample + "_R1.fastq.gz"

    elif os.path.isfile(pairedR1) and read == "reverse":
        return sample_dir + sample + "_R2.fastq.gz"

    #elif os.path.isfile(paired1) and read == "forward":
    #    return sample_dir + sample + "_1.fastq.gz"

    #elif os.path.isfile(paired1) and read == "reverse":
    #    return sample_dir + sample + "_2.fastq.gz"
    
    else:
        return ""


def input4filter(wildcards, sampleList, srrEnabled, srrList):

    if isSingle(wildcards.samples, sampleList, srrEnabled, srrList):
        return "results/{samples}/{samples}_{build}_se.bed"
    else:    
        return "results/{samples}/{samples}_{build}_pe.bed"

def input4PCA(sampleList, srrEnabled, srrList, build):

    inputList = []
    for sample in sampleList:

        if isSingle(sample, sampleList, srrEnabled, srrList):
            inputList.append("results/" + sample + "/" + sample + "_" + build + "_se_sortedbyCoordinates.bam")
        else:    
            inputList.append("results/" + sample + "/" + sample + "_" + build + "_pe_sortedbyCoordinates.bam")

    return inputList

def input4nucTable(method):

    if method == "XR":
        return "results/{samples}/{samples}_{build}_lengthMode.fa"
    elif method == "DS":    
        return "results/{samples}/{samples}_{build}_sorted_10.fa"

def getDamage(sample, damage_type, sampleList):

    try:
        idx = sampleList.index(sample)
    except:
       raise(ValueError("Designated wildcard cannot be found in sample list."))
        
    return damage_type[idx]

def getMotif(sample, damageList, sampleList):
    
    tDamage = getDamage(sample, damageList, sampleList)

    if tDamage.lower() in ["oxaliplatin", "cisplatin", "bpdedg"]: 
        return "'.{4}(g|G){2}.{4}'"
    
    elif tDamage.lower() in ["64", "64pp", "(6-4)pp", "6-4pp", "cpd"]: 
        return "'.{4}(c|t|C|T){2}.{4}'"

    elif tDamage.lower() == "na":
        return "'.{10}'"

def getDinuc(sample, damageList, sampleList):

    tDamage = getDamage(sample, damageList, sampleList)
    
    if tDamage.lower() in ["oxaliplatin", "cisplatin"]: 
        return "'GG'"
    
    elif tDamage.lower() in ["64", "64pp", "(6-4)pp", "6-4pp", "cpd"]: 
        return "'CC','CT','TC','TT'"

def getInput(sample, inputExist, inputList, inputIdx, sampleList):

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
            
                return k
                
    else:
        return ""

def lineNum(file):
    
    linenum = 0
    if os.path.exists(file):
        with open(file) as f:
            for line in f:
                linenum += 1

    warnMessage = ("\n" + file + " file is either empty or does not exists!\n" + 
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

def allInput(method, build, sampleList, srrEnabled, srrList):

    inputList = []
    for sample in sampleList:
        sampledir = "results/" + sample + "/" 

        if isSingle(sample, sampleList, srrEnabled, srrList):
            inputList.append(sampledir + sample + ".html")
        else:
            inputList.append(sampledir + sample + "_1.html")
            inputList.append(sampledir + sample + "_2.html")

        inputList.append(sampledir + sample + "_" + build + 
            "_sorted_nucleotideTable.png")
        inputList.append(sampledir + sample + "_" + build + 
            "_sorted_dinucleotideTable.png")
        inputList.append(sampledir + sample + "_" + build + "_" + method +
            "_sim_nucleotideTable.png")
        inputList.append(sampledir + sample + "_" + build + "_" + method + 
            "_sim_dinucleotideTable.png")
        inputList.append(sampledir + sample + "_" + build + 
            "_length_distribution.png")
        inputList.append(sampledir + sample + "_" + build + "_" + method + 
            "_sorted_plus.bw")
        inputList.append(sampledir + sample + "_" + build + "_" + method + 
            "_sorted_minus.bw")
        #inputList.append(sampledir + sample + "_" + build + 
        #    "_igv_report_chrX.html")

        #for i in range(1,23):
        #    inputList.append(sampledir + sample + "_" + build + 
        #        "_igv_report_chr" + str(i) + ".html")

        if method == "DS":
            inputList.append(sampledir + sample + "_" + build + 
                "_sorted_ds_dipyrimidines_plus.bed") 
            inputList.append(sampledir + sample + "_" + build + 
                "_sorted_ds_dipyrimidines_minus.bed") 

        elif method == "XR":
            inputList.append(sampledir + sample + "_" + 
                build + "_sorted_plus.bed") 
            inputList.append(sampledir + sample + "_" + 
                build + "_sorted_minus.bed") 
    
    if method == "DS":
        inputList.append("results/scatterplot_PearsonCorr_bigwigScores_DS.png")
        inputList.append("results/PearsonCorr_bigwigScores_DS.tab")
        inputList.append("results/heatmap_SpearmanCorr_readCounts_DS.png")
        inputList.append("results/SpearmanCorr_readCounts_DS.tab")
        inputList.append("results/PCA_readCounts_DS.png")
    elif method == "XR":
        inputList.append("results/scatterplot_PearsonCorr_bigwigScores_XR.png")
        inputList.append("results/PearsonCorr_bigwigScores_XR.tab")
        inputList.append("results/heatmap_SpearmanCorr_readCounts_XR.png")
        inputList.append("results/SpearmanCorr_readCounts_XR.tab")
        inputList.append("results/PCA_readCounts_XR.png")
        
    #print(inputList)
    return inputList

################################################################################
