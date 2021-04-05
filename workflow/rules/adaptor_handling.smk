
rule cutadapt_se:
    input:
        "resources/samples/{samples}.fastq.gz",
    output:
        fastq="results/{samples}/{samples}_cutadapt.fastq.gz",
        qc="results/{samples}/{samples}_cutadapt.qc.txt",   
    params:
        adapters=config["adaptor_se"],  
        extra=config["cutadapt_se"],  
    log:
        "logs/{samples}/{samples}_cutadapt.log",
    benchmark:
        "logs/{samples}/{samples}_cutadapt.benchmark.txt",
    wrapper:
        "0.69.0/bio/cutadapt/se"

rule cutadapt_pe:
    input:
        fq1="resources/samples/{samples}_R1.fastq.gz",
        fq2="resources/samples/{samples}_R2.fastq.gz",
    output:
        fastq1="results/{samples}/{samples}_cutadapt_1.fastq.gz",
        fastq2="results/{samples}/{samples}_cutadapt_2.fastq.gz",
        qc="results/{samples}/{samples}_cutadapt.qc.txt",
    params:
        adapters=config["adaptor_pe"],  
        extra=config["cutadapt_pe"], 
    log:
        "logs/{samples}/{samples}_cutadapt.log",
    benchmark:
        "logs/{samples}/{samples}_cutadapt.benchmark.txt",
    wrapper:
        "0.69.0/bio/cutadapt/pe"