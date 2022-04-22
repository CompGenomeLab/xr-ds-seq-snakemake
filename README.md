# XR-seq and Damage-seq workflows

This repository contains xr-seq and damage-seq workflows.  

<br>

## Installation

- This workflow is prepared using 
    [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management 
    system and [conda](https://docs.conda.io/en/latest/)

- To run the workflow, you should have conda installed for environment 
    management. All the other packages including Snakemake and their 
    dependencies can be obtained automatically through environments prepared 
    for each step of the workflow. You can follow the installation steps from 
    [the link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html).

- Initially, you should clone the repository and navigate into the directory: 

    ```
    git clone https://github.com/CompGenomeLab/xr-ds-seq-snakemake.git
        
    cd xr-ds-seq-snakemake
    ```

- Next, you should create a conda environment with the defined packages. 
    Install [mamba](https://mamba.readthedocs.io/en/latest/) 
    and create the environment using mamba:

    ```
    conda install -c conda-forge mamba

    mamba create -c bioconda -c conda-forge -c r -n repair snakemake=6.3.0 python=3.8 rust=1.55.0

    conda activate repair
    ```

<br>

## Install Simulation

- You should install [boquila](https://github.com/CompGenomeLab/boquila) 
    tool to simulate repair and damage reads.

- Path of the program should be added to the PATH variable:

    ```
    export PATH="$PATH:~/.cargo/bin" # works on linux
    ```

<br>

## Directory Structure

This workflow is prepared according to the 
[structure](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html) 
recommended by Snakemake: 

- `config/`: contains the configuration files.

- `logs/`: contains the log files of each step. 
    This folder will automatically appear when you run the workflow.

- `report/`: contains the description files of figures,
    which will be used in reports.

- `resources/`: contains `samples/` where the raw XR-seq and Damage-seq data 
    are stored, `input/` where the input files are stored, 
    and `ref_genomes/` where the reference genome files are stored. 
    Reference genome files can be automatically produced by the workflows, 
    if they are properly defined in the config files.  

- `results/`: contains the generated files and figures. 
    This folder will automatically appear when you run the workflow.

- `workflow/`: contains `envs/` where the environments are stored,
    `rules/` where the Snakemake rules are stored, and 
    `scripts/` where the scripts used inside the rules are stored. 

<br>

## Configuration file

Before running the workflow, you should edit the configuration files. 
For both XR-seq and Damage-seq workflows, there are 2 configuration files: 
`config_(XR/DS)_initial.yaml` and `config_(XR/DS).yaml`. 
In most of the cases, 
the configuration file with "initial" prefix shouldn't be modified 
by the user since they are containing configuration settings 
that are common for all XR-seq and Damage-seq experiments. 
A config example and description of each parameter 
for "config_(XR/DS).yaml" are given below:

```
sample: 
 - "NHF1_CPD_1h_XR_rep1"
 - "NHF1_CPD_1h_XR_rep2"

meta: 
  NHF1_CPD_1h_XR_rep1:
    srr_id: "SRR3062593:SRR3062594:SRR3062595" 
    layout: "single"
    product: "CPD"
    simulation_enabled: True
    simulation_input: "SRR5461463" 
    simulation_input_layout: "paired"
  NHF1_CPD_1h_XR_rep2:
    srr_id: "SRR3062596:SRR3062597:SRR3062598" 
    layout: "single"
    product: "CPD"
    simulation_enabled: True
    simulation_input: "SRR5461463" 
    simulation_input_layout: "paired"

genome:
  build: "hg19"
  link: "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz"
```

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
        `{SAMPLE}_R1.fastq.gz`, `{SAMPLE}_R2.fastq.gz` and 
        `{SAMPLE}_1.fastq.gz`, `{SAMPLE}_2.fastq.gz` as raw data.
        Therefore, paired-end sample files must contain `_R1/2` or `_1/2` 
        suffixes and the suffixes should not be given 
        in the configuration file to the `sample`.

- `meta`: contains the metadata of each sample. 

    - `srr_id`: The SRR code of the sample. 

        - Each downloaded raw data with the SRR codes will be named as the 
            corresponding sample name in the `sample` parameter.  

        - If a sample have multiple SRR codes, then the code should be provided 
            in the given format: `SRRXXXXXXX:SRRXXXXXXX:SRRXXXXXXX` 

        - If the fastq file is already provided in the `resources/samples/` 
            directory, workflow will directly use that file. In such a case,
            you don't have to provide this parameter.

    - `layout`: Whether the given sample is sequenced as 
        paired-end or single-end. 
        'Single' or 'Paired' (Case-insensitive)

    - `product`: Damage type of each sample. Currently damages below are 
        available can be provided as (case-insensitive):

        - (6-4)PP: `64`, `64pp`, `(6-4)pp`, `6-4pp`;
        - CPD: `CPD`;
        - Cisplatin: `cisplatin`;
        - Oxaliplatin: `oxaliplatin`.

    <br>    

    - `simulation_enabled`: True if you want to generate simulated reads of 
        your sample via [boquila](https://github.com/CompGenomeLab/boquila). 
        If not, False should be provided.

    - `simulation_input`: Simulation can be either done using 
        the reference genome or with a provided input file. If an input file
        will be used, it's SRR code should be given. 

        - If the srr file of the input is already provided in the 
            `resources/input/` directory in gzipped fastq format, 
            workflow will directly use that file. Even in that case,
            this parameter should be used as it will be used to set the name
            of the input file.

        - If simulation will be done without an input file, this parameter 
            should be removed.

    - `simulation_input_layout`: Whether the given input file is sequenced as 
        paired-end or single-end. 
        'Single' or 'Paired' (Case-insensitive)

        - If simulation will be done without an input file, this parameter 
            should be removed.

- `genome`: Contains 2 parameters, which are `build` and `link`.

    - `build`: The name of the reference genome that the workflow will use.
        Any name desired by the user can be used.

    - `link`: The url of the reference fasta file to be retrieved. The file 
        should be gzipped.

<br>

## Usage

After adjusting the configuration file, you can run the workflow 
from `xr-ds-seq-snakemake` directory.

- For XR-seq:

    ```
    snakemake -pr --use-conda --cores 64 --snakefile Snakefile_XR.smk --keep-going --rerun-incomplete 
    ```

- For Damage-seq:

    ```
    snakemake -pr --use-conda --cores 64 --snakefile Snakefile_DS.smk --keep-going --rerun-incomplete 
    ```
| Note: To run the workflow on [Slurm Workload Manager](https://slurm.schedmd.com/srun.html) as set of jobs, `--profile` flag must be provided: |  
| --- |
    snakemake -pr --use-conda --profile config/slurm --snakefile Snakefile_DS.smk --keep-going --rerun-incomplete 

<br>

To generate detailed HTML report files, 
the code below should be run after workflow:

```
snakemake --report report.zip
```