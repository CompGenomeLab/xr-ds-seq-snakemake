
rule plot_length:
    input:
        rules.length_distribution.output
    output:
        report("results/processed_files/{samples}_{build}_{method}_length_distribution.pdf", 
                category="Length Distribution"),
    params:
        "{samples}_{build}_{method}",
    log:
        "logs/rule/plot_length/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/plot_length/{samples}_{build}_{method}.benchmark.txt",
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