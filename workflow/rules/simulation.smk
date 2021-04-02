
rule simulation_ds:
    input:
        plus="results/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        minus="results/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed",
        genome="resources/ref_genomes/{build}/genome_{build}.fa",
        regions="resources/ref_genomes/{build}/genome_{build}.ron", 
    output:
        bed=temp("results/{samples}/{samples}_{build}_sorted_ds_dipyrimidines.bed"),
        fa=temp("results/{samples}/{samples}_{build}_sorted_ds_dipyrimidines.fa"),
        sim="results/{samples}/{samples}_{build}_ds_sim.fa",     
    log:
        "logs/{samples}/{samples}_{build}_simulation_ds.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_simulation_ds.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        cat {input.plus} {input.minus} > {output.bed}
        
        bedtools getfasta \
        -fi {input.genome} \
        -bed {output.bed} \
        -fo {output.fa} \
        -s

        boquila \
        --fasta {output.fa} \
        --ref {input.genome} \
        --regions {input.regions} \
        > {output.sim}
        """

rule simulation_xr:
    input:
        bed="results/{samples}/{samples}_{build}_sorted_chr.bed",
        genome="resources/ref_genomes/{build}/genome_{build}.fa",
        regions="resources/ref_genomes/{build}/genome_{build}.ron", 
    output:
        fa=temp("results/{samples}/{samples}_{build}_sorted_chr.fa"),
        sim="results/{samples}/{samples}_{build}_xr_sim.fa",
    log:
        "logs/{samples}/{samples}_{build}_simulation_xr.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_simulation_xr.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.bed} \
        -fo {output.fa} \
        -s

        boquila \
        --fasta {output.fa} \
        --ref {input.genome} \
        --regions {input.regions} \
        > {output.sim}
        """