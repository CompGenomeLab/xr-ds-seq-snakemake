
import warnings
from os import path

################### Helper Functions ###########################################

def isSingle(sample):

    sample_dir = "resources/samples/"
    single = sample_dir + sample + ".fastq.gz"
    pairedR1 = sample_dir + sample + "_R1.fastq.gz"
    
    if path.isfile(single) and "_R1" not in single:
        return True
    elif path.isfile(pairedR1):
        return False
    else:
        raise(ValueError("Sample not found..."))

def getMotif(wildcards):
    
    if "oxaliplatin" in wildcards.samples or "cisplatin" in wildcards.samples: 
        return "'.{4}(g|G){2}.{4}'"
    
    elif "64" in wildcards.samples or "CPD" in wildcards.samples:
        return "'.{4}(c|t|C|T){2}.{4}'"

def getDinuc(wildcards):
    
    if "Oxaliplatin" in wildcards.samples or "Cisplatin" in wildcards.samples: 
        return "'GG'"
    
    elif "64" in wildcards.samples or "CPD" in wildcards.samples or "R190" in wildcards.samples:
        return "'CC','CT','TC','TT'"

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

def allInput(method, build, sampleList):

    inputList = []
    for sample in sampleList:
        sampledir = "results/" + sample + "/" 
    
        if isSingle(sample):
            inputList.append(sampledir + sample + ".html")
        else:
            inputList.append(sampledir + sample + "_R1.html")
            inputList.append(sampledir + sample + "_R2.html")

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
                
    return inputList

################################################################################