
rule sort_filter:
    input:
        "results/{dir}/{samples}{v}/{samples}_cutadapt.bed",
    output:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_chr.bed",
    params:
        filt=lambda w: getSampleInfo(w, config["filter"]),
    log:
        "results/{dir}/{samples}{v}/log/sort_filter.log",
    benchmark:
        "results/{dir}/{samples}{v}/log/sort_filter.benchmark.txt",
    shell:  
        """
        sort -u -k1,1 -k2,2n -k3,3n {input} |
        egrep {params.filt} > {output}
        """


