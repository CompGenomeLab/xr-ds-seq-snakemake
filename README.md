# XR-seq and Damage-seq pipelines

This repository contains xr-seq and damage-seq pipelines.  

## Installation

- This workflow is prepared using 
[Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management 
system and [conda](https://docs.conda.io/en/latest/)

- To run the workflow, you should have conda installed for environment 
management. All the other packages including Snakemake and their dependencies 
can be obtained automatically through environments prepared for each step of 
the workflow. You can follow the installation steps from 
[the link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html).

- Initially, you should clone the repository and navigate into the directory: 

    ```
    git clone https://github.com/CompGenomeLab/xr-ds-seq-snakemake.git
        
    cd xr-ds-seq-snakemake
    ```

- Next, you should create a conda environment with the defined packages. 
We propose 2 way to create the environment:

    - One is installing [mamba](https://mamba.readthedocs.io/en/latest/) 
    and creating the environment using mamba:

        ```
        conda install -c conda-forge mamba

        mamba create -c bioconda -c conda-forge -c r -n repair snakemake python=3.8 rust=1.50

        conda activate repair
        ```

    - Or the environment can be directly created from our environment file:

        ```
        conda env create -f workflow/envs/env.yaml

        conda activate repair
        ```

## Install Simulation

- You should install [boquila](https://github.com/CompGenomeLab/boquila) 
tool to simulate repair and damage reads.

    ```
    cargo install --branch main --git https://github.com/CompGenomeLab/boquila.git
    ```

| Warning: please add the path of the tool to your $PATH to be able to run boquila tool. |
| --- |

## Directory Structure (will be updated)

This workflow is prepared according to the 
[structure](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html) 
recommended by Snakemake: 

- `config/`: contains the configuration file.
- `resources/`: contains the input files.  
- `results/`: contains the generated files and figures. 
This folder will automatically appear when you run the workflow.
- `logs/`: contains the log files of each step. 
This folder will automatically appear when you run the workflow.
- `workflow/`: contains `envs/` where the environments are stored, 
`rules/` where the Snakemake rules are stored, and 
`scripts/` where the scripts used inside the rules are stored. 

## Configuration file (will be updated)

Before running the workflow, you should edit the configuration file.  

## Usage (will be updated)

After adjusting the configuration file, you can run the workflow 
from `xr-ds-seq-snakemake` directory:

```
snakemake -pr --use-conda --cores 64 --snakefile Snakefile_(XR/DS) --debug-dag
```
