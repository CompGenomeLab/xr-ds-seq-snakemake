
rule plot_nuc_sim:
    input:
        dinuc="results/{method}/{samples}/{samples}_{build}_{method}_sim_dinucleotideTable.txt",
        nuc="results/{method}/{samples}/{samples}_{build}_{method}_sim_nucleotideTable.txt",
    output:
        dinuc=report("results/{method}/{samples}/{samples}_{build}_{method}_sim_dinucleotideTable.pdf", 
                    category="Nucleotide Content"),
        nuc=report("results/{method}/{samples}/{samples}_{build}_{method}_sim_nucleotideTable.pdf", 
                    category="Nucleotide Content"),
    params:
        motif=lambda w: getDinuc(w.samples, config["meta"][w.samples]["product"]),
        name="{samples}_{build}_{method}_sim",
    log:
        "logs/rule/figs/{samples}/{samples}_{build}_{method}_plot_nuc_sim.log",
    benchmark:
        "logs/rule/figs/{samples}/{samples}_{build}_{method}_plot_nuc_sim.benchmark.txt",
    conda:
        "../envs/plot_nuc.yaml"
    shell:  
        """
        Rscript workflow/scripts/nucleotidePlot.r \
        -i {input.dinuc} \
        -k 2 \
        -s {params.name} \
        -f {params.motif} \
        -o {output.dinuc} \
        -l {log}

        Rscript workflow/scripts/nucleotidePlot.r \
        -i {input.nuc} \
        -k 1 \
        -s {params.name} \
        -o {output.nuc} \
        -l {log}
        """