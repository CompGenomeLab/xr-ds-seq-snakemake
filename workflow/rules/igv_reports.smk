
rule igv_report:
    input:
        fasta="resources/ref_genomes/{build}/genome_{build}_chr{chr}.fa",
        vcf="resources/ref_genomes/{build}/variation_{build}_bgzip.vcf.gz",
        tracks=["results/{samples}/{samples}_{build}_sorted_plus.bed", "results/{samples}/{samples}_{build}_sorted_minus.bed"],
    output:
        report("results/{samples}/{samples}_{build}_igv_report_chr{chr}.html", category="IGV"),
    params:
        extra="",  
    log:
        "logs/{samples}/{samples}_{build}_igv_report_chr{chr}.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_igv_report_chr{chr}.benchmark.txt",
    wrapper:
        "0.69.0/bio/igv-reports"