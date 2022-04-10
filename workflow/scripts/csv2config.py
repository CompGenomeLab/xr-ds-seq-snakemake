
import csv
import json
import os

def defineMethod(methodList):

    if ("Damage-seq" or "DS") in methodList:
        return "DS"
    elif ("XR-seq" or "XR") in methodList:
        return "XR"

def removeAfterLastRegex(mystring, myregex):

    if myregex in mystring:
        return (mystring[::-1].split(myregex[::-1])[1])[::-1]

def checkVersionNum(configList):
    
    method = defineMethod(configList[1])
    dirList = []
    version = 1

    for root, dirs, files in os.walk("results/" + method):
        for mydir in dirs:
            if mydir.startswith(configList[2]):
                dirList.append(mydir)
    #print(dirList)

    for mydir in dirList:
        filename = "results/" + method + "/" + mydir + "/version.info.txt" 

        with open(filename, "r") as f:
            fread = f.readline()
            if fread.split("', '")[3:] == str(configList).split("', '")[3:]:
                return fread.split("', '")[2].split(".v")[1]
        
        version += 1
        #print(version)
    
    return str(version)

def createVersionFile(configList):

    method = defineMethod(configList[1])
    filename = "results/" + method + "/" + configList[2] + "/version.info.txt"
    os.makedirs(os.path.dirname(filename), exist_ok=True)

    with open(filename, "w") as f:
        f.write(str(configList))

def updateVersion(configList):

    method = defineMethod(configList[1])
    checkVersion = False

    for root, dirs, files in os.walk("results/" + method):
        for folder in dirs:
            folder = removeAfterLastRegex(folder, ".v")
            if configList[2] == folder:
                checkVersion = True
                break
    
    if checkVersion:   
        versionNum = checkVersionNum(configList)
        configList[2] += (".v" + versionNum)
        createVersionFile(configList)
        return configList

    else:
        configList[2] += ".v1"
        createVersionFile(configList)
        return configList

def updateDict(listInfo, headerList, configDict):
    
    for idx in range(2,len(listInfo)):
        if headerList[idx] in configDict:
            configDict[headerList[idx]].append(listInfo[idx])
        else:
            configDict[headerList[idx]] = [listInfo[idx]]
    
    return configDict

def duplicateSamples(listInfo):

    sampleInfoList = [[]]
    for info in listInfo:
        semicolon = info.count(";")
        infoList = [info]

        if semicolon > 0:
            infoList = info.strip().split(";")
        
            tempList = []
            for spltInfo in infoList:
                for subList in sampleInfoList:
                    subListemp = subList[:]
                    subListemp.append(spltInfo)
                    tempList.append(subListemp)

            sampleInfoList = [] 
            for subList in tempList:
                sampleInfoList.append(subList)

        else:
            for subList in sampleInfoList:
                subList.append(info)

    return sampleInfoList

def csv2Dict(csvfile, methodList):

    configDict = {}
    configDict["method"] = defineMethod(methodList)

    with open(csvfile) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:

            if row[0] == "no":
                headerList = row

            if row[1] in methodList:
                dupRows = duplicateSamples(row)
                for dupRow in dupRows:
                    dupRow = updateVersion(dupRow)
                    configDict = updateDict(dupRow, headerList, configDict)
        
    return configDict

def dict2json(config, myjson):

    with open(myjson, "w") as outfile:  
        json.dump(config, outfile)

#methodList = ["DS", "Damage-seq"]

#config = csv2Dict("config/samples.csv", methodList)

#print(checkSamples("HDA64A1_ATCACG_sub10k", "DS"))

#print(removeAfterLastRegex("hello.vhadi.vneolacak", ".v"))

#myList = ['1', 'Damage-seq', 'HDA_ATCACG_sub10k', 'single-end', '-g GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT', '--discard-trimmed', '', 'GRCh37', '-q 20', "'^([1-9]|1[0-9]|2[0-2]|X)'", '.{4}(c|t|C|T){2}.{4}']

#checkVersionNum(myList)

#print(removeAfterLastRegex("cemo.v23asd", ".v"))