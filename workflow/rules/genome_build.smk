
rule genome_build:
    input:
        reference="resources/ref_genomes/{build}/genome_{build}.fa",
    output:
        multiext(
        "resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    log: 
        "resources/ref_genomes/{build}/log/bowtie2_build.log"
    benchmark:
        "resources/ref_genomes/{build}/log/bowtie2_build.benchmark.txt",
    params:
        extra=""  
    threads: 
        4
    wrapper:
        "0.69.0/bio/bowtie2/build"