rule mark_duplicates_se:
    input:
        "results/{method}/{samples}/{samples}_cutadapt_se_{build}_samSorted.bam"
    output:
        bam="results/{method}/{samples}/{samples}_dedup_cutadapt_se_{build}.bam",
        metrics="results/{method}/{samples}/{samples}_cutadapt_se_dedup_{build}.metrics.txt",
    log:
        "logs/rule/analysis/{samples}/{samples}_{method}_{build}_se_mark_duplicates.log",
    params:
        extra="--REMOVE_DUPLICATES true",
    wrapper:
        "v1.3.2/bio/picard/markduplicates"

rule mark_duplicates_pe:
    input:
        "results/{method}/{samples}/{samples}_cutadapt_pe_{build}_samSorted.bam"
    output:
        bam="results/{method}/{samples}/{samples}_dedup_cutadapt_pe_{build}.bam",
        metrics="results/{method}/{samples}/{samples}_cutadapt_pe_dedup_{build}.metrics.txt",
    log:
        "logs/rule/analysis/{samples}/{samples}_{method}_{build}_pe_mark_duplicates.log",
    params:
        extra="--REMOVE_DUPLICATES true",
    wrapper:
        "v1.3.2/bio/picard/markduplicates"