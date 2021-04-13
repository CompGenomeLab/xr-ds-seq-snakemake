
rule fastqc:
    input:
        "resources/samples/{samples}.fastq.gz",
    output:
        html=report("results/{samples}/{samples}.html", category="QC"),
        zip="results/{samples}/{samples}_fastqc.zip",
    params: ""
    log:
        "logs/{samples}/{samples}_fastqc.log",
    benchmark:
        "logs/{samples}/{samples}_fastqc.benchmark.txt",
    threads: 1
    wrapper:
        "0.69.0/bio/fastqc"