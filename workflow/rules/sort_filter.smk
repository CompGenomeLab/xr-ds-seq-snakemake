
rule sort_filter:
    input:
        "results/{samples}/{samples}_{build}_.bed",
    output:
        "results/{samples}/{samples}_{build}_sorted_chr.bed",
    params:
        filt=config["filter"],
    log:
        "logs/{samples}/{samples}_{build}_sort_filter.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_sort_filter.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Sorting and filtering bed file by chromosomes..." &&
        sort -u -k1,1 -k2,2n -k3,3n {input} |&
        egrep {params.filt} > {output} &&
        echo "`date -R`: Success! Bed file is filtered." || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        """