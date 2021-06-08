
rule genome_chr_download:
    output:
        "resources/ref_genomes/{build}/genome_{build}_chr{chr}.fa"
    params:
        species=config["species"],
        datatype=config["datatype"],
        build=config["build"],
        release=config["release"],
        chromosome="{chr}",
    log:
        "resources/ref_genomes/{build}/log/download_chr{chr}.log",
    benchmark:
        "resources/ref_genomes/{build}/log/download_chr{chr}.benchmark.txt",
    wrapper:
        "0.69.0/bio/reference/ensembl-sequence"