
rule bed2fasta_ds:
    input:
        plus="results/{samples}/{samples}_{build}_sorted_plus_10.bed",
        minus="results/{samples}/{samples}_{build}_sorted_minus_10.bed",
        genome="resources/ref_genomes/{build}/genome_{build}.fa",
    output:
        plus="results/{samples}/{samples}_{build}_sorted_plus_10.fa",
        minus="results/{samples}/{samples}_{build}_sorted_minus_10.fa",
        comb="results/{samples}/{samples}_{build}_sorted_10.fa",
        bed=temp("results/{samples}/{samples}_{build}_sorted_10.bed"),       
    log:
        "logs/{samples}/{samples}_{build}_bed2fasta_ds.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_bed2fasta_ds.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        cat {input.plus} {input.minus} > {output.bed}
        
        bedtools getfasta -fi {input.genome} -bed {output.bed} \
        -fo {output.comb} -s

        bedtools getfasta -fi {input.genome} -bed {input.plus} \
        -fo {output.plus} -s

        bedtools getfasta -fi {input.genome} -bed {input.minus} \
        -fo {output.minus} -s
        """

rule bed2fasta_xr:
    input:
        bed="results/{samples}/{samples}_{build}_lengthMode.bed",
        genome="resources/ref_genomes/{build}/genome_{build}.fa",
    output:
        "results/{samples}/{samples}_{build}_lengthMode.fa",
    log:
        "logs/{samples}/{samples}_{build}_bed2fasta_xr.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_bed2fasta_xr.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        bedtools getfasta -fi {input.genome} -bed {input.bed} -fo {output} -s
        """