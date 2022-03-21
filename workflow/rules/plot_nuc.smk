
rule plot_nuc:
    input:
        dinuc="results/{method}/{samples}/{samples}_{build}_sorted_dinucleotideTable.txt",
        nuc="results/{method}/{samples}/{samples}_{build}_sorted_nucleotideTable.txt",
    output:
        dinuc=report("results/{method}/{samples}/{samples}_{build}_sorted_dinucleotideTable.pdf", 
                    category="Nucleotide Content"),
        nuc=report("results/{method}/{samples}/{samples}_{build}_sorted_nucleotideTable.pdf", 
                    category="Nucleotide Content"),
    params:
        motif=lambda w: getDinuc(w.samples, config["damage_type"], config["sample"]),
        name="{samples}_{build}",
    log:
        "logs/rule/figs/{samples}/{samples}_{build}_{method}_plot_nuc.log",
    benchmark:
        "logs/rule/figs/{samples}/{samples}_{build}_{method}_plot_nuc.benchmark.txt",
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