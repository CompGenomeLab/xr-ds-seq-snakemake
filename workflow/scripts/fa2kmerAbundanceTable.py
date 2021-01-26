#!/usr/bin/env python

import fasta

if snakemake.params[2]:
    percentageFlag = True
else:
    percentageFlag = False
    
Fasta = fasta.fasta(snakemake.input[0])
if snakemake.params[1]:
    kmerAbundanceDict = Fasta.getKmerAbundance(int(snakemake.params[0]), snakemake.params[1])
else:
    kmerAbundanceDict = Fasta.getKmerAbundance(int(snakemake.params[0]))

Fasta.writeKmerAbundanceTable(kmerAbundanceDict, snakemake.output[0], percentageFlag)
