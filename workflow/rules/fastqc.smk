
rule fastqc_se:
    input:
        "resources/samples/{samples}.fastq.gz",
    output:
        html=report("results/processed_files/{samples}.html", category="QC"),
        zip="results/processed_files/{samples}_fastqc.zip",
    wildcard_constraints:
        samples='|'.join([s for s in config["meta"].keys()]),
    params: ""
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_fastqc.benchmark.txt",
    threads: 1
    wrapper:
        "0.69.0/bio/fastqc"

rule fastqc_pe:
    input:
        "resources/samples/{samples}_{ext}.fastq.gz", 
    output:
        html=report("results/processed_files/{samples}_{ext}.html", category="QC"), 
        zip="results/processed_files/{samples}_{ext}_fastqc.zip", 
    wildcard_constraints:
        samples='|'.join([s for s in config["meta"].keys()]),
    params: ""
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_fastqc_{ext}.benchmark.txt",
    threads: 1
    wrapper:
        "0.69.0/bio/fastqc"