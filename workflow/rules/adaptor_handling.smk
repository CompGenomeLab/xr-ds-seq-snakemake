
rule cutadapt_se:
    input:
        rules.sra_se.output,
    output:
        fastq=temp("results/{method}/{samples}/{samples}_cutadapt.fastq.gz"),
        qc=report("results/{method}/{samples}/{samples}_cutadapt.qc.txt", category="QC"),   
    params:
        adapters=lambda w: getMethodParams(w, config["meta"], "adaptor", 
            config["XR"], config["DS"]),
        extra=lambda w: getMethodParams(w, config["meta"], "cutadapt", 
            config["XR"], config["DS"]),
    log:
        "logs/rule/analysis/{samples}/{samples}_{method}_cutadapt.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{method}_cutadapt.benchmark.txt",
    wrapper:
        "0.75.0/bio/cutadapt/se"

rule cutadapt_pe:
    input:
        fq1=rules.sra_pe.output.read1, 
        fq2=rules.sra_pe.output.read2,
    output:
        fastq1=temp("results/{method}/{samples}/{samples}_cutadapt_1.fastq.gz"),
        fastq2=temp("results/{method}/{samples}/{samples}_cutadapt_2.fastq.gz"),
        qc=report("results/{method}/{samples}/{samples}_cutadapt.qc.txt", category="QC"),
    params:
        adapters=lambda w: getMethodParams(w, config["meta"], "adaptor", 
            config["XR"], config["DS"]),
        extra=lambda w: getMethodParams(w, config["meta"], "cutadapt", 
            config["XR"], config["DS"]),
    log:
        "logs/rule/analysis/{samples}/{samples}_{method}_cutadapt.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{method}_cutadapt.benchmark.txt",
    wrapper:
        "0.75.0/bio/cutadapt/pe"
