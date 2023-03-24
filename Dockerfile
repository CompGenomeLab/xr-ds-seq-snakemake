FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="2795c90a388a0867ec21c476d74184b5ada7f0b0c9badeaf86b71756586b9166"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/align.yaml
#   prefix: /conda-envs/43837407976759d77bf1da9da340da71
#   name: align
#   dependencies:
#     - bioconda::bowtie2 ==2.4.1  
#     - bioconda::samtools ==1.10
#     - conda-forge::tbb ==2020.2
RUN mkdir -p /conda-envs/43837407976759d77bf1da9da340da71
COPY workflow/envs/align.yaml /conda-envs/43837407976759d77bf1da9da340da71/environment.yaml

# Conda environment:
#   source: workflow/envs/bam_correlation.yaml
#   prefix: /conda-envs/17ca4a013d053833603718e62172e570
#   name: bam_correlation
#   dependencies:
#     - bioconda::deeptools ==3.1.3
RUN mkdir -p /conda-envs/17ca4a013d053833603718e62172e570
COPY workflow/envs/bam_correlation.yaml /conda-envs/17ca4a013d053833603718e62172e570/environment.yaml

# Conda environment:
#   source: workflow/envs/bedGraphToBigWig.yaml
#   prefix: /conda-envs/b3dc2e81133eede24d82db429706a643
#   name: bedGraphToBigWig
#   dependencies:
#     - bioconda::ucsc-bedgraphtobigwig == 377
RUN mkdir -p /conda-envs/b3dc2e81133eede24d82db429706a643
COPY workflow/envs/bedGraphToBigWig.yaml /conda-envs/b3dc2e81133eede24d82db429706a643/environment.yaml

# Conda environment:
#   source: workflow/envs/bedtools.yaml
#   prefix: /conda-envs/c897f951bd717f91534019efe1296dd6
#   name: bedtools
#   dependencies:
#     - bioconda::samtools ==1.10
#     - bioconda::bedtools =2.29.0
RUN mkdir -p /conda-envs/c897f951bd717f91534019efe1296dd6
COPY workflow/envs/bedtools.yaml /conda-envs/c897f951bd717f91534019efe1296dd6/environment.yaml

# Conda environment:
#   source: workflow/envs/cutadapt.yaml
#   prefix: /conda-envs/904dc933e84239ba20c1956456252a84
#   name: cutadapt
#   dependencies:
#     - bioconda::cutadapt =4.1
RUN mkdir -p /conda-envs/904dc933e84239ba20c1956456252a84
COPY workflow/envs/cutadapt.yaml /conda-envs/904dc933e84239ba20c1956456252a84/environment.yaml

# Conda environment:
#   source: workflow/envs/fastqc.yaml
#   prefix: /conda-envs/626ff9f2bedbfc6368271be101cfe6fc
#   name: fastqc
#   dependencies:
#     - bioconda::fastqc ==0.11.9
RUN mkdir -p /conda-envs/626ff9f2bedbfc6368271be101cfe6fc
COPY workflow/envs/fastqc.yaml /conda-envs/626ff9f2bedbfc6368271be101cfe6fc/environment.yaml

# Conda environment:
#   source: workflow/envs/picard.yaml
#   prefix: /conda-envs/4b4f4dd9bad1debd8d48c405b4a6415e
#   name: picard
#   dependencies:
#     - bioconda::picard =2.27
#     - bioconda::samtools =1.15
RUN mkdir -p /conda-envs/4b4f4dd9bad1debd8d48c405b4a6415e
COPY workflow/envs/picard.yaml /conda-envs/4b4f4dd9bad1debd8d48c405b4a6415e/environment.yaml

# Conda environment:
#   source: workflow/envs/plot_nuc.yaml
#   prefix: /conda-envs/095ac112bf1947f3cf84f022c319c8fb
#   name: plot_nuc
#   dependencies:
#     - conda-forge::r-base ==3.5.1
#     - conda-forge::r-argparser ==0.4
#     - conda-forge::r-ggplot2 ==3.0.0
#     - conda-forge::r-dplyr ==0.8.0.1
#     - conda-forge::r-tidyr ==0.8.2
#     - conda-forge::r-futile.logger ==1.4.3
RUN mkdir -p /conda-envs/095ac112bf1947f3cf84f022c319c8fb
COPY workflow/envs/plot_nuc.yaml /conda-envs/095ac112bf1947f3cf84f022c319c8fb/environment.yaml

# Conda environment:
#   source: workflow/envs/simulation.yaml
#   prefix: /conda-envs/b95925071e90477716cb31b521bb84a9
#   name: simulation
#   dependencies:
#     - bioconda::bedtools =2.29.0
#     - bioconda::bowtie2 ==2.4.1  
#     - bioconda::samtools ==1.10
#     - conda-forge::tbb ==2020.2
#     - bioconda::boquila ==0.6.0
RUN mkdir -p /conda-envs/b95925071e90477716cb31b521bb84a9
COPY workflow/envs/simulation.yaml /conda-envs/b95925071e90477716cb31b521bb84a9/environment.yaml

# Conda environment:
#   source: workflow/envs/sra.yaml
#   prefix: /conda-envs/25ef06dcfd9a3944b25499c543688a44
#   name: sra
#   dependencies:
#     - bioconda::sra-tools ==3.0.0
RUN mkdir -p /conda-envs/25ef06dcfd9a3944b25499c543688a44
COPY workflow/envs/sra.yaml /conda-envs/25ef06dcfd9a3944b25499c543688a44/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/43837407976759d77bf1da9da340da71 --file /conda-envs/43837407976759d77bf1da9da340da71/environment.yaml && \
    mamba env create --prefix /conda-envs/17ca4a013d053833603718e62172e570 --file /conda-envs/17ca4a013d053833603718e62172e570/environment.yaml && \
    mamba env create --prefix /conda-envs/b3dc2e81133eede24d82db429706a643 --file /conda-envs/b3dc2e81133eede24d82db429706a643/environment.yaml && \
    mamba env create --prefix /conda-envs/c897f951bd717f91534019efe1296dd6 --file /conda-envs/c897f951bd717f91534019efe1296dd6/environment.yaml && \
    mamba env create --prefix /conda-envs/904dc933e84239ba20c1956456252a84 --file /conda-envs/904dc933e84239ba20c1956456252a84/environment.yaml && \
    mamba env create --prefix /conda-envs/626ff9f2bedbfc6368271be101cfe6fc --file /conda-envs/626ff9f2bedbfc6368271be101cfe6fc/environment.yaml && \
    mamba env create --prefix /conda-envs/4b4f4dd9bad1debd8d48c405b4a6415e --file /conda-envs/4b4f4dd9bad1debd8d48c405b4a6415e/environment.yaml && \
    mamba env create --prefix /conda-envs/095ac112bf1947f3cf84f022c319c8fb --file /conda-envs/095ac112bf1947f3cf84f022c319c8fb/environment.yaml && \
    mamba env create --prefix /conda-envs/b95925071e90477716cb31b521bb84a9 --file /conda-envs/b95925071e90477716cb31b521bb84a9/environment.yaml && \
    mamba env create --prefix /conda-envs/25ef06dcfd9a3944b25499c543688a44 --file /conda-envs/25ef06dcfd9a3944b25499c543688a44/environment.yaml && \
    mamba clean --all -y
