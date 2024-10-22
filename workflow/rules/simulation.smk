
rule simulation_ds:
    input:
        bed=rules.comb_strands.output,
        genome=rules.genome_download.output,
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2",
        inp_file=lambda w: input4rule(w, config["meta"], "simulation_ds"),
    output:
        bed="results/{method}/{samples}/{samples}_{build}_DS_sim.bed",  
        fa="results/{method}/{samples}/{samples}_{build}_ds_sim.fa",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
    threads: 16
    log:
        "logs/rule/simulation_ds/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/simulation_ds/{samples}_{build}_{method}.benchmark.txt",
    conda:
        "../envs/simulation.yaml"
    shell:
        """
        workflow/scripts/simulation.sh \
        {input.bed} \
        {input.genome} \
        {input.inp_file} \
        {params.ref_genome} \
        {output.bed} \
        {output.fa} \
        {log}
        """

rule simulation_xr:
    input:
        bed=rules.sort_filter.output,
        genome=rules.genome_download.output,
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2", 
        inp_file=lambda w: input4rule(w, config["meta"], "simulation_xr"),
    output:
        bed="results/{method}/{samples}/{samples}_{build}_XR_sim.bed",
        fa="results/{method}/{samples}/{samples}_{build}_xr_sim.fa",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
    threads: 16
    log:
        "logs/rule/simulation_xr/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/simulation_xr/{samples}_{build}_{method}.benchmark.txt",
    conda:
        "../envs/simulation.yaml"
    shell:
        """
        workflow/scripts/simulation.sh \
        {input.bed} \
        {input.genome} \
        {input.inp_file} \
        {params.ref_genome} \
        {output.bed} \
        {output.fa} \
        {log}
        """