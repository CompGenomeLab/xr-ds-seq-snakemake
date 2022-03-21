
rule genome_index2ron:
    input:
        "resources/ref_genomes/{build}/genome_{build}.fa.fai",
    output:
        "resources/ref_genomes/{build}/genome_{build}.ron",
    log:
        "logs/rule/analysis/{build}/log/indexing.log",
    benchmark:
        "logs/rule/analysis/{build}/log/indexing.benchmark.txt",
    shell:
        "python3 workflow/scripts/idx2ron.py -i {input} -o {output} -l {log}"