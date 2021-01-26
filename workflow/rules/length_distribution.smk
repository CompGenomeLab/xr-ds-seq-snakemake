
rule length_dist:
    input:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_chr.bed"
    output:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_length_distribution.txt"
    log:
        "results/{dir}/{samples}{v}/log/length_dist.log"
    benchmark:
        "results/{dir}/{samples}{v}/log/length_dist.benchmark.txt",
    shell:  
        """
        awk '{{print $3-$2}}' {input} | sort -k1,1n | uniq -c | 
        sed 's/\s\s*/ /g' | awk '{{print $2"\t"$1}}' > {output}
        """









