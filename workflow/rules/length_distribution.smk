
rule length_dist:
    input:
        "results/{samples}/{samples}_{build}_sorted_chr.bed",
    output:
        "results/{samples}/{samples}_{build}_length_distribution.txt",
    log:
        "logs/{samples}/{samples}_{build}_length_dist.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_length_dist.benchmark.txt",
    shell:  
        """
        awk '{{print $3-$2}}' {input} | sort -k1,1n | uniq -c | 
        sed 's/\s\s*/ /g' | awk '{{print $2"\\t"$1}}' > {output}
        """









