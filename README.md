# XR-seq and Damage-seq workflows

This repository contains xr-seq and damage-seq workflows.  

<br>

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
Install [mamba](https://mamba.readthedocs.io/en/latest/) 
and create the environment using mamba:

    ```
    conda install -c conda-forge mamba

    mamba create -c bioconda -c conda-forge -c r -n repair snakemake=6.3.0 python=3.8 rust=1.50 sra-tools=2.11.0

    conda activate repair
    ```

<br>

## Install Simulation

- You should install [boquila](https://github.com/CompGenomeLab/boquila) 
tool to simulate repair and damage reads.

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
    `{SAMPLE}_R1.fastq.gz`, `{SAMPLE}_R2.fastq.gz` and 
    `{SAMPLE}_1.fastq.gz`, `{SAMPLE}_2.fastq.gz` as raw data.
    Therefore, paired-end sample files must contain `_R1/2` or `_1/2` suffixes and 
    the suffixes should not be given in the configuration file to the `sample`.

- `damage_type`: Damage type of each sample should be provided here in the 
same order of the samples:

    ```
    damage_type: 
        - "64"
        - "CPD"
        - "oxaliplatin"
    ```

    - Currently damages below are available can be provided as (case-insensitive):

        - (6-4)PP: `64`, `64pp`, `(6-4)pp`, `6-4pp`;
        - CPD: `CPD`;
        - Cisplatin: `cisplatin`;
        - Oxaliplatin: `oxaliplatin`.

- `srr`: Contains 2 parameters: `enabled`, `codes`. `enabled` is a boolean 
variable (can only be True or False) and if it is True, then the pipeline will
retrieve the raw data by [sra-toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc) 
the defining the SRR code in the `codes` parameter. 

    ```
    srr: 
    enabled: True
    codes:
        - "SRRXXXXXXX"
    ```

    - Each downloaded raw data with the SRR codes will be named as the 
    corresponding sample name in the `sample` parameter.  

    - If a sample have multiple SRR codes, then the code should be provided in 
    the given format: `SRRXXXXXXX:SRRXXXXXXX:SRRXXXXXXX` 

- `input`: Contains 4 parameters: `exist`, `files`, `sample`, `srr`. This parameter is used for [boquila](https://github.com/CompGenomeLab/boquila), which can take input files for more accurate read simulations. If you don't have any input file, you can set `exist` as False. `files` parameter contains the names of input files w/o the extension. `sample` parameter contains the indeces of samples that will use the input of the same order. Lastly, `srr` can be used for downloading inputs from [sra-toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc). It works the same was as the `srr` parameter above. As a whole inputs can be provided as:

    ```
    input:
        exist: True
        files:
            - "input0"
            - "input1"
            - "input2"  
        sample: 
            - "0-3" # input0 will be used for sample 0,1,2,3
            - "4, 7" # input1 will be used for sample 4,7
            - "5, 6" # input2 will be used for sample 5,6
        srr: 
            enabled: True
            codes:
            - "SRRXXXXXXX" # srr of input0
            - "SRRXXXXXXX" # srr of input1
            - "SRRXXXXXXX" # srr of input2
    ```

- Genome parameters: The parameters `build`, `species`, `datatype`, 
and `release` are defined to retrieve correct reference genome from ensembl. 
For more information, you can check out the 
[link](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/reference/ensembl-sequence.html). 

    - If `genome_download` is set to False, then a fasta file expected to be 
    provided as `resources/ref_genomes/${build}/genome_${build}.fa`, where 
    `${build}` should be the string given in `build` parameter. 
    
    - In the same manner, If `bowtie2_build` is set to False, 
    then a build should be provided in `resources/ref_genomes/${build}/Bowtie2/`, 
    where `${build}` should be the string given in `build` parameter. 
    File name before the extensions must be `genome_${build}`. Lastly to 
    properly set `bowtie2_build` to False, `genome_download` must be set 
    to False as well.

- `filter`: This parameter is used to filter chromosomes by the given regex.

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