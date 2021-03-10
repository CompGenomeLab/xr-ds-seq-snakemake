#!/usr/bin/env python

import fasta
import argparse
import re
import sys

parser = argparse.ArgumentParser(description='converts fasta to bed by choosing a motif of interest')
parser.add_argument('-i', required= False, help='input')
parser.add_argument('-o', required= False, help='output')
parser.add_argument('-r', required= True, help='regex motif that is expected')

args = parser.parse_args()

fastaObject = fasta.fasta(args.i)
output = args.o
regexMotif = args.r
pattern = re.compile(regexMotif)

def criterionPass(citerion):
    code = 'True if ' + citerion + ' else False' 
    result = eval(code)
    return result

def fastaHeader2bedLine(header):
    chr = header.replace('>','').split(':')[0]
    start = int(header.split(':')[1].split('(')[0].split('-')[0])
    end = int(header.split(':')[1].split('(')[0].split('-')[1])
    strand = header.split(':')[1].split('(')[1][0]
    return chr + '\t' + str(start) + '\t' + str(end) + '\t.\t.\t' + strand

out = open(output, 'w')
seqDicts = fastaObject.stream(100*4096)
for seqDict in seqDicts:
    header = seqDict['h']
    sequence = seqDict['s']
    if pattern.match(sequence):
        out.write(fastaHeader2bedLine(header) + '\n')
out.close()