
import os
import os.path
import warnings
from os import path

################### Helper Functions ###########################################

def getSampleInfo(wildcards, info):
    mysample = wildcards.samples + wildcards.v
    idx = config["sample"].index(mysample)
    return info[idx]
    
def getGenomeInfo(wildcards, info):
    idx = config["build"].index(wildcards.build)
    return info[idx]

def getGenome(wildcards, gtype="Bowtie2"):
    mysample = wildcards.samples + wildcards.v
    idx = config["sample"].index(mysample)
    build=config["build"][idx]

    if gtype == "Bowtie2":
        return "resources/ref_genomes/" + build + "/Bowtie2/genome_" + build
    elif gtype == "genome":
        return "resources/ref_genomes/" + build + "/genome_" + build + ".fa"
    elif gtype == "index":    
        return "resources/ref_genomes/" + build + "/genome_" + build + ".fa.fai"

def checkRef(build):
    
    ref_genome = "resources/ref_genomes/" + build + "/Bowtie2/genome_" + build
    extensions = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", 
    ".rev.2.bt2"]
    ref_genome_list = []

    for ext in extensions:
        ref_genome_list.append(ref_genome + ext)
    
    for genome in ref_genome_list:
        if not path.exists(genome):
            return False

    return True

def RefAlert(build):

    missingBuilds = ""
    for genome in build:
        if checkRef(genome) == False:
            if genome not in missingBuilds:
                missingBuilds += genome + ", "

    if len(missingBuilds) > 0:
        missingBuilds = missingBuilds[:-2]
        errMessage = ("{0} build/s is/are missing or has/have missing files. " + 
            "Please download {0} build/s properly before the pipeline " + 
            "execution.").format(missingBuilds)

        raise ValueError(errMessage)

def lineNum(file):
    
    linenum = 0
    if path.exists(file):
        with open(file) as f:
            for line in f:
                linenum += 1

    warnMessage = ("\n" + file + " file is either empty or does not exists!\n" + 
        "It is expected if this is a dry-run. The file will be produced " + 
        "after the execution.")

    if linenum == 0:
        warnings.warn(warnMessage)

    return linenum

def mappedReads(file):

    file1 = str(file[0])
    file2 = str(file[1])

    return lineNum(file1) + lineNum(file2)

def allInput(method, sampleList):

    inputList = []
    for sampleV in sampleList:
        sampledir = "results/" + method + "/" + sampleV + "/"

        if ".v" in sampleV:
            sample = (sampleV[::-1].split(".v"[::-1])[1])[::-1]

        if method == "DS":
            inputList.append(sampledir + sample + 
                "_cutadapt_sorted_plus_dipyrimidines.bw")
            inputList.append(sampledir + sample + 
                "_cutadapt_sorted_minus_dipyrimidines.bw")
            inputList.append(sampledir + sample + 
                "_cutadapt_sorted_dinucleotideTable.txt")
            inputList.append(sampledir + sample + 
                "_cutadapt_length_distribution.txt")

        elif method == "XR": 
            inputList.append(sampledir + sample + "_cutadapt_sorted_plus.bw")
            inputList.append(sampledir + sample + "_cutadapt_sorted_minus.bw")
            inputList.append(sampledir + sample + 
                "_cutadapt_sorted_dinucleotideTable.txt")
            inputList.append(sampledir + sample + 
                "_cutadapt_length_distribution.txt")
    
    return inputList



################################################################################
