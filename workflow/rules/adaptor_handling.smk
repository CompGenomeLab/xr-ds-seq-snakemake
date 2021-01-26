
rule cutadapt_se:
    input:
        "resources/samples/{samples}.fastq.gz",
    output:
        fastq="results/{dir}/{samples}{v}/{samples}_cutadapt.fastq.gz",
        qc="results/{dir}/{samples}{v}/log/cutadapt.qc.txt",   
    params:
        adapters=lambda w: getSampleInfo(w, config["adaptor"]),  
        extra=lambda w: getSampleInfo(w, config["cutadapt"]),  
    log:
        "results/{dir}/{samples}{v}/log/cutadapt.log",
    benchmark:
        "results/{dir}/{samples}{v}/log/cutadapt.benchmark.txt",
    wrapper:
        "0.69.0/bio/cutadapt/se"

rule cutadapt_pe:
    input:
        fq1="resources/samples/{samples}_R1.fastq.gz",
        fq2="resources/samples/{samples}_R2.fastq.gz",
    output:
        fastq1="results/{dir}/{samples}{v}/{samples}_cutadapt_1.fastq.gz",
        fastq2="results/{dir}/{samples}{v}/{samples}_cutadapt_2.fastq.gz",
        qc="results/{dir}/{samples}{v}/log/cutadapt.qc.txt",
    params:
        adapters=lambda w: getSampleInfo(w, config["adaptor"]),  
        extra=lambda w: getSampleInfo(w, config["cutadapt"]), 
    log:
        "results/{dir}/{samples}{v}/log/cutadapt.log",
    benchmark:
        "results/{dir}/{samples}{v}/log/cutadapt.benchmark.txt",
    wrapper:
        "0.69.0/bio/cutadapt/pe"
