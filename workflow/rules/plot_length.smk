
rule plot_length:
    input:
        "results/{method}/{samples}/{samples}_{build}_length_distribution.txt",
    output:
        report("results/{method}/{samples}/{samples}_{build}_length_distribution.pdf", 
                category="Length Distribution"),
    params:
        "{samples}_{build}",
    log:
        "logs/rule/figs/{samples}/{samples}_{build}_{method}_plot_length.log",
    benchmark:
        "logs/rule/figs/{samples}/{samples}_{build}_{method}_plot_length.benchmark.txt",
    conda:
        "../envs/plot_nuc.yaml"
    shell:  
        """
        Rscript workflow/scripts/lengthDistPlot.r \
        -i {input} \
        -s {params} \
        -o {output} \
        -l {log}
        """