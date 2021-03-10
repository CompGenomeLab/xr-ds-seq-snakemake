
rule filtbyMotifs:
    input:
        "results/{samples}/{samples}_{build}_sorted_{strand}_10.fa",
    output:
        "results/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_{strand}.bed",
    params:
        lambda wc: config["motif"],
    log:
        "logs/{samples}/{samples}_{build}_filtbyMotifs_{strand}.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_filtbyMotifs_{strand}.benchmark.txt",
    shell:  
        """
        workflow/scripts/fa2bedByChoosingReadMotifs.py \
        -i {input} \
        -o {output} \
        -r {params}
        """

