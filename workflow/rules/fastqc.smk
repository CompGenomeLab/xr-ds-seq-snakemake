
rule fastqc_se:
    input:
        "resources/samples/{samples}.fastq.gz",
    output:
        html=report("results/{method}/{samples}/{samples}.html", category="QC"),
        zip="results/{method}/{samples}/{samples}_fastqc.zip",
    params: ""
    log:
        "logs/rule/analysis/{samples}/{samples}_{method}_fastqc.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{method}_fastqc.benchmark.txt",
    threads: 1
    wrapper:
        "0.69.0/bio/fastqc"

rule fastqc_pe:
    input:
        "resources/samples/{samples}_{ext}.fastq.gz", 
    output:
        html=report("results/{method}/{samples}/{samples}_{ext}.html", category="QC"), 
        zip="results/{method}/{samples}/{samples}_{ext}_fastqc.zip", 
    params: ""
    log:
        "logs/rule/analysis/{samples}/{samples}_{method}_fastqc_{ext}.log", 
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{method}_fastqc_{ext}.benchmark.txt",
    threads: 1
    wrapper:
        "0.69.0/bio/fastqc"