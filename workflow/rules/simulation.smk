
rule simulation_ds:
    input:
        bed=rules.comb_strands.output,
        genome=rules.genome_download.output,
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2",
        inp_file=lambda w: getInput(w.samples, config["meta"], config["genome"]["build"]),
    output:
        bed="results/{method}/{samples}/{samples}_{build}_ds_sim.bed",  
        fa="results/{method}/{samples}/{samples}_{build}_ds_sim.fa",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_simulation_ds.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_simulation_ds.benchmark.txt",
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
        inp_file=lambda w: getInput(w.samples, config["meta"], config["genome"]["build"]),
    output:
        bed="results/{method}/{samples}/{samples}_{build}_xr_sim.bed",
        fa="results/{method}/{samples}/{samples}_{build}_xr_sim.fa",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_simulation_xr.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_simulation_xr.benchmark.txt",
    conda:
        "../envs/simulation.yaml"
    shell:
        """
        workflow/scripts/simulation.sh \
        {input.bed} \
        {input.genome} \
        {input.inp_file} \
        {params.ref_genome} \
        {output.fa} \
        {output.bed} \
        {log}
        """