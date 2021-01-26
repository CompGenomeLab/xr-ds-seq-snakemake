
rule bedGraphToBigWig_ds:
    input:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_{strand}_dipyrimidines.bdg",
    output:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_{strand}_dipyrimidines.bw",
    params:
        index=lambda w: getGenome(w, "index"),
    log:
        "results/{dir}/{samples}{v}/log/bedGraphToBigWig_{strand}.log",
    benchmark:
        "results/{dir}/{samples}{v}/log/bedGraphToBigWig_{strand}.benchmark.txt",
    conda:
        "../envs/bedGraphToBigWig.yaml"
    shell:  
        """
        bedGraphToBigWig {input} {params.index} {output}
        """

rule bedGraphToBigWig_xr:
    input:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_{strand}.bdg",
    output:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_{strand}.bw",
    params:
        index=lambda w: getGenome(w, "index"),
    log:
        "results/{dir}/{samples}{v}/log/bedGraphToBigWig_{strand}.log",
    benchmark:
        "results/{dir}/{samples}{v}/log/bedGraphToBigWig_{strand}.benchmark.txt",
    conda:
        "../envs/bedGraphToBigWig.yaml"
    shell:  
        """
        bedGraphToBigWig {input} {params.index} {output}
        """