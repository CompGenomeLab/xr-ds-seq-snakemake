
rule plot_nuc_sim:
    input:
        dinuc="results/{samples}/{samples}_{build}_sim_dinucleotideTable.txt",
        nuc="results/{samples}/{samples}_{build}_sim_nucleotideTable.txt",
    output:
        dinuc=report("results/{samples}/{samples}_{build}_sim_dinucleotideTable.png"),
        nuc=report("results/{samples}/{samples}_{build}_sim_nucleotideTable.png"),
    params:
        motif=lambda w: getDinuc(w),
        name="{samples}_{build}_sim",
    log:
        "logs/{samples}/{samples}_{build}_plot_nuc_sim.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_plot_nuc_sim.benchmark.txt",
    conda:
        "../envs/plot_nuc.yaml"
    shell:  
        """
        Rscript workflow/scripts/nucleotidePlot.r \
        -i {input.dinuc} \
        -k 2 \
        -s {params.name} \
        -f {params.motif} \
        -o {output.dinuc} 2>{log}

        Rscript workflow/scripts/nucleotidePlot.r \
        -i {input.nuc} \
        -k 1 \
        -s {params.name} \
        -o {output.nuc} 2>>{log}
        """