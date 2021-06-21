
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
        pairedR1 = sample_dir + sample + "_1.fastq.gz"
        
        if os.path.isfile(single) and "_1" not in single:
            return True
        elif os.path.isfile(pairedR1):
            return False
        else:
            raise(ValueError("Sample not found..."))

def input4filter(wildcards, sampleList, srrEnabled, srrList):

    if isSingle(wildcards.samples, sampleList, srrEnabled, srrList):
        return "results/{samples}/{samples}_{build}_se.bed"
    else:    
        return "results/{samples}/{samples}_{build}_pe.bed"

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

    if tDamage.lower() in ["oxaliplatin", "cisplatin"]: 
        return "'.{4}(g|G){2}.{4}'"
    
    elif tDamage.lower() in ["64", "64pp", "(6-4)pp", "6-4pp", "cpd"]: 
        return "'.{4}(c|t|C|T){2}.{4}'"

def getDinuc(wildcards):
    
    if tDamage.lower() in ["oxaliplatin", "cisplatin"]: 
        return "'GG'"
    
    elif tDamage.lower() in ["64", "64pp", "(6-4)pp", "6-4pp", "cpd"]: 
        return "'CC','CT','TC','TT'"

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

def mappedReads(file):

    file1 = str(file[0])
    file2 = str(file[1])

    return lineNum(file1) + lineNum(file2)

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
    
    #print(inputList)
    return inputList

################################################################################