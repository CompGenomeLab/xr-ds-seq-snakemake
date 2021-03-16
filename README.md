# xr-ds-seq-snakemake
This repository contains xr-seq and damage-seq pipelines.  



conda install -c conda-forge mamba

mamba create -c bioconda -c conda-forge -c r -n repair snakemake python=3.8 rust=1.50

conda activate repair

## add simulation

cargo install --branch main --git https://github.com/CompGenomeLab/boquila.git

Warning: please add the path of the tool to your $PATH to be able to run boquila tool.

