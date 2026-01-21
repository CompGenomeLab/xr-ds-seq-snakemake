#!/bin/bash

#SBATCH -J combinelanes

#SBATCH --account=adelab
#SBATCH --qos=adelab
#SBATCH --partition=genomics

#SBATCH --output=/cta/users/cazgari/pipelines/xr-ds-seq-snakemake/logs/cluster/combinelanes-%j.out
#SBATCH --error=/cta/users/cazgari/pipelines/xr-ds-seq-snakemake/logs/cluster/combinelanes-%j.err

#SBATCH --array=0-40

output_path="/cta/users/cazgari/pipelines/xr-ds-seq-snakemake/resources/samples"
sample_path="/cta/users/cazgari/pipelines/xr-ds-seq-snakemake/resources/samples"

filelist=($(find $sample_path -maxdepth 1 -type f -name "*L001_R1_001.fastq.gz" | awk -F '/' '{print $NF}' | sed 's/_L001_R1_001\.fastq\.gz$//' | sort -u))
sample=${filelist[$SLURM_ARRAY_TASK_ID]}

echo "Processing sample: $sample"
echo "Sample slurm array ID: $SLURM_ARRAY_TASK_ID"

zcat $sample_path/${sample}_L001_R1_001.fastq.gz $sample_path/${sample}_L002_R1_001.fastq.gz > $output_path/${sample}_R1.fastq

gzip $output_path/${sample}_R1.fastq

rm $sample_path/${sample}_L001_R1_001.fastq.gz
rm $sample_path/${sample}_L002_R1_001.fastq.gz