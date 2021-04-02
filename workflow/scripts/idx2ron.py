#!/usr/bin/env python

import argparse
import re
import sys

parser = argparse.ArgumentParser(description='converts genome index file to ron file.')
parser.add_argument('-i', required= False, help='input')
parser.add_argument('-o', required= False, help='output')

args = parser.parse_args()

idx_file = args.i
output = args.o

with open(idx_file, "r") as f:
    lines = f.readlines()

chrList = []
for line in lines:
    ll = line.strip().split("\t")
    newline = '(name: "' + ll[0] + '", start: 0, end: ' + ll[1] + ')'
    chrList.append(newline)

with open(output, 'w') as fout:
    fout.write("[\n")
    for item in chrList:
        fout.write("%s,\n" % item)
    fout.write("]")