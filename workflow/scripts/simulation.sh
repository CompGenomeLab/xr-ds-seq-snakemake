#!/usr/bin/env bash

set -e 

(echo "`date -R`: Converting {input.bed} to fasta format..." &&
bedtools getfasta \
-fi $2 \
-bed $1 \
-fo $1.fa \
-s &&
echo "`date -R`: Success! {input.bed} is converted." || 
{{ echo "`date -R`: Process failed..."; exit 1; }}  ) > $7 2>&1

if [ "$3" == *.ron ]; then

    (echo "`date -R`: Simulating reads according to reference genome..." &&
    boquila \
    --fasta $1.fa \
    --bed $5 \
    --ref $2 \
    --seed 1 \
    --sens 2 \
    --regions $3 \
    > $6 &&
    {{ echo "`date -R`: Success! Simulation is done."; rm $1.fa }} || 
    {{ echo "`date -R`: Process failed..."; rm $5; rm $6; exit 1; }} \
    ) >> $7 2>&1

else

    (echo "`date -R`: Simulating reads with input file..." &&
    boquila \
    --fasta $1.fa \
    --inseqFasta \
    --inseq $3 \
    --seed 1 \
    --sens 2 \
    > $6 &&
    {{ echo "`date -R`: Success! Simulation is done."; rm $1.fa }}  || 
    {{ echo "`date -R`: Process failed..."; rm $6; exit 1; }}  ) >> $7 2>&1

    (echo "`date -R`: Aligning fasta file..." &&
    bowtie2 \
    -x $4 \
    -f $6 -S $6.sam &&
    echo "`date -R`: Success! Alignment is done." || 
    {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> $7 2>&1

    (echo "`date -R`: Converting sam to bam..." &&
    samtools view -Sbh -o $6.bam $6.sam &&
    {{ echo "`date -R`: Success! Conversion is done."; rm $6.sam }} || 
    {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> $7 2>&1

    (echo "`date -R`: Sorting (coordinates) bam file..." &&
    samtools sort $6.bam > $6_sorted.bam &&
    {{ echo "`date -R`: Success! Bam file is sorted."; rm $6.bam }} || 
    {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> $7 2>&1

    (echo "`date -R`: Index bam file..." &&
    samtools index $6_sorted.bam $6_sorted.bam.bai &&
    echo "`date -R`: Success! Bam file is sorted." || 
    {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> $7 2>&1

    (echo "`date -R`: Processing bam file..." && 
    samtools view -b $6_sorted.bam |&
    bedtools bamtobed > $5 &&
    {{echo "`date -R`: Success! Bam file converted to bed format."; \
    rm $6_sorted.bam; rm $6_sorted.bam.bai }} || 
    {{ echo "`date -R`: Process failed..."; rm $5; exit 1; }}  ) >> $7 2>&1

fi