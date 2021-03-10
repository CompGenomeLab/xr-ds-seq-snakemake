
rule genome_indexing:
    input:
        "resources/ref_genomes/{build}/genome_{build}.fa",
    output:
        "resources/ref_genomes/{build}/genome_{build}.fa.fai",
    benchmark:
        "resources/ref_genomes/{build}/log/indexing.benchmark.txt",
    wrapper:
        "0.69.0/bio/samtools/faidx"