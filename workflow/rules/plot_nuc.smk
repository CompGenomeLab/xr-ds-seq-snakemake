
rule plot_nuc:
    input:
        dinuc="results/{samples}/{samples}_{build}_sorted_dinucleotideTable.txt",
        nuc="results/{samples}/{samples}_{build}_sorted_nucleotideTable.txt",
    output:
        dinuc=report("results/{samples}/{samples}_{build}_sorted_dinucleotideTable.png", 
                    category="Nucleotide Content"),
        nuc=report("results/{samples}/{samples}_{build}_sorted_nucleotideTable.png", 
                    category="Nucleotide Content"),
    params:
        motif=lambda w: getDinuc(w),
        name="{samples}_{build}",
    log:
        "logs/{samples}/{samples}_{build}_plot_nuc.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_plot_nuc.benchmark.txt",
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