
rule igv_report:
    input:
        fasta="resources/ref_genomes/{build}/genome_{build}.fa",
        vcf="resources/ref_genomes/{build}/variation_{build}_bgzip.vcf.gz",
        tracks=["results/{samples}/{samples}_{build}_sorted_plus.bed", "results/{samples}/{samples}_{build}_sorted_minus.bed"],
    output:
        "results/{samples}/{samples}_{build}_igv_report.html",
    params:
        extra="",  
    log:
        "logs/{samples}/{samples}_{build}_igv_report.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_igv_report.benchmark.txt",
    wrapper:
        "0.69.0/bio/igv-reports"