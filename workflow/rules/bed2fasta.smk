
rule bed2fasta_ds:
    input:
        plus="results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_plus_10.bed",
        minus="results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_minus_10.bed",
    output:
        plus="results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_plus_10.fa",
        minus="results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_minus_10.fa",
        comb="results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_10.fa",
        bed=temp("results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_10.bed"),
    params:
        genome=lambda w: getGenome(w, "genome"),      
    log:
        "results/{dir}/{samples}{v}/log/bed2fasta_ds.log",
    benchmark:
        "results/{dir}/{samples}{v}/log/bed2fasta_ds.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        cat {input.plus} {input.minus} > {output.bed}
        
        bedtools getfasta -fi {params.genome} -bed {output.bed} \
        -fo {output.comb} -s

        bedtools getfasta -fi {params.genome} -bed {input.plus} \
        -fo {output.plus} -s

        bedtools getfasta -fi {params.genome} -bed {input.minus} \
        -fo {output.minus} -s
        """

rule bed2fasta_xr:
    input:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_lengthMode.bed",
    output:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_lengthMode.fa",
    params:
        genome=lambda w: getGenome(w, "genome"),       
    log:
        "results/{dir}/{samples}{v}/log/bed2fasta_xr.log",
    benchmark:
        "results/{dir}/{samples}{v}/log/bed2fasta_xr.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        bedtools getfasta -fi {params.genome} -bed {input} -fo {output} -s
        """