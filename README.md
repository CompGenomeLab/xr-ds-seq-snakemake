# XR-seq and Damage-seq workflows

This repository contains xr-seq and damage-seq workflows.  

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

        mamba create -c bioconda -c conda-forge -c r -n repair snakemake python=3.8 rust=1.50 sra-tools=2.11.0

        conda activate repair
        ```

    - Or the environment can be directly created from our environment file:

        ```
        conda env create -f workflow/envs/env.yaml

        conda activate repair
        ```

## Install Simulation

- You should install [boquila](https://github.com/CompGenomeLab/boquila) 
tool to simulate repair and damage reads:

    ```
    cargo install --branch main --git https://github.com/CompGenomeLab/boquila.git
    ```

| Warning: please add the path of the tool to your $PATH to be able to run boquila tool. |
| --- |

## Directory Structure

This workflow is prepared according to the 
[structure](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html) 
recommended by Snakemake: 

- `config/`: contains the configuration files.

- `logs/`: contains the log files of each step. 
This folder will automatically appear when you run the workflow.

- `reports/`: contains the report files, which can be produced 
after the workflow is over. 

- `resources/`: contains `samples/` where the raw XR-seq and Damage-seq data 
are stored and `ref_genomes/` where the reference genome files are stored. 
Reference genome files can be automatically produced by the workflows, 
if they are properly defined in the config files.  

- `results/`: contains the generated files and figures. 
This folder will automatically appear when you run the workflow.

- `workflow/`: contains `envs/` where the environments are stored, 
`rules/` where the Snakemake rules are stored, and 
`scripts/` where the scripts used inside the rules are stored. 

## Configuration file

Before running the workflow, you should edit the configuration files. 
For both XR-seq and Damage-seq workflows, there are 2 configuration files: 
`config_(XR/DS)_initial.yaml` and `config_(XR/DS).yaml`. 
The configuration file with "_initial_" prefix shouldn't be modified 
by the user since they are containing configuration settings 
that are common for all XR-seq and Damage-seq experiments. 
For more detail about these configuration files, 
check out the readme file in `config/` directory. 
The parameters for "config_(XR/DS).yaml" as below:

- `sample`: The name of the sample file w/o the extension. 
Multiple sample names can be given in the below format:

    ```
    sample: 
    - "SAMPLE_1"
    - "SAMPLE_2"
    - "SAMPLE_3"
    ```

    - Using the given sample name, the workflow will look for 
    `{SAMPLE}.fastq.gz` as raw data. 
    Therefore, the fastq file must be gzipped before running the workflow.

    - If the layout of the given sample is paired-end, 
    the workflow will look for 
    `{SAMPLE}_R1.fastq.gz` and `{SAMPLE}_R2.fastq.gz` as raw data.
    Therefore, paired-end sample files must contain `_R1/2` suffixes and 
    the suffixes should not be given in the configuration file to the `sample`.

    - Lastly, because damage type is retrieved from the name of the file, 
    given sample name should contain the damage type. 
    If damage type is:

        - (6-4)PP, then sample name should contain `64`.
        - CPD, then sample name should contain `CPD`.
        - Cisplatin, then sample name should contain `cisplatin`.
        - Oxaliplatin, then sample name should contain `oxaliplatin`.

- Genome parameters: The parameters `build`, `species`, `datatype`, 
and `release` are defined to retrieve correct reference genome from ensembl. 
For more information, you can check out the 
[link](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/reference/ensembl-sequence.html). 

- `filter`: This parameter is used to filter chromosomes by the given regex.

## Usage

After adjusting the configuration file, you can run the workflow 
from `xr-ds-seq-snakemake` directory.

- For XR-seq:

    ```
    snakemake -pr --use-conda --cores 64 --snakefile Snakefile_XR --debug-dag
    ```

- For Damage-seq:

    ```
    snakemake -pr --use-conda --cores 64 --snakefile Snakefile_DS --debug-dag
    ```

To generate detailed HTML report files, 
the code below should be run after workflow:

```
snakemake -pr --use-conda --cores 64 --snakefile Snakefile_(XR/DS) --debug-dag --report reports/{report_name}.html
```