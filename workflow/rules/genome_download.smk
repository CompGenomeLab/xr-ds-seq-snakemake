
rule genome_download:
    output:
        "resources/ref_genomes/{build}/genome_{build}.fa"
    params:
        species=config["species"],
        datatype=config["datatype"],
        build=config["build"],
        release=config["release"],
    log:
        "resources/ref_genomes/{build}/log/download.log",
    benchmark:
        "resources/ref_genomes/{build}/log/download.benchmark.txt",
    wrapper:
        "0.69.0/bio/reference/ensembl-sequence"
