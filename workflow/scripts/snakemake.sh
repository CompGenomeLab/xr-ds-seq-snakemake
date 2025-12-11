#!/bin/bash

#SBATCH -J snakemake

#SBATCH --account=adelab
#SBATCH --qos=adelab
#SBATCH --partition=genomics

#SBATCH --error=/cta/users/cazgari/pipelines/xr-ds-seq-snakemake/logs/cluster/snakemake-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

eval "$(conda shell.bash hook)"
conda activate snakemake

snakemake --unlock
snakemake --profile config/slurm
