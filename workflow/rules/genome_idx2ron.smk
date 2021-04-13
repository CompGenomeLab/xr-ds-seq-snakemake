
rule genome_index2ron:
    input:
        "resources/ref_genomes/{build}/genome_{build}.fa.fai",
    output:
        "resources/ref_genomes/{build}/genome_{build}.ron",
    log:
        "resources/ref_genomes/{build}/log/indexing.log",
    benchmark:
        "resources/ref_genomes/{build}/log/indexing.benchmark.txt",
    shell:
        "python3 workflow/scripts/idx2ron.py -i {input} -o {output} -l {log}"