
rule filtbyMotifs:
    input:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_{strand}_10.fa",
    output:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_{strand}_dipyrimidines.bed",
    params:
        motif=lambda w: getSampleInfo(w, config["motif"])
    log:
        "results/{dir}/{samples}{v}/log/filtbyMotifs_{strand}.log",
    benchmark:
        "results/{dir}/{samples}{v}/log/filtbyMotifs_{strand}.benchmark.txt",
    script:  
        "../scripts/fa2bedByChoosingReadMotifs.py"

