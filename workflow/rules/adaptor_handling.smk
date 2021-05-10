
rule cutadapt_se:
    input:
        "resources/samples/{samples}.fastq.gz",
    output:
        fastq=temp("results/{samples}/{samples}_cutadapt.fastq.gz"),
        qc=report("results/{samples}/{samples}_cutadapt.qc.txt", category="QC"),   
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
        fq1="resources/samples/{samples}_1.fastq.gz",
        fq2="resources/samples/{samples}_2.fastq.gz",
    output:
        fastq1=temp("results/{samples}/{samples}_cutadapt_1.fastq.gz"),
        fastq2=temp("results/{samples}/{samples}_cutadapt_2.fastq.gz"),
        qc=report("results/{samples}/{samples}_cutadapt.qc.txt", category="QC"),
    params:
        adapters=config["adaptor_pe"],  
        extra=config["cutadapt_pe"], 
    log:
        "logs/{samples}/{samples}_cutadapt.log",
    benchmark:
        "logs/{samples}/{samples}_cutadapt.benchmark.txt",
    wrapper:
        "0.69.0/bio/cutadapt/pe"