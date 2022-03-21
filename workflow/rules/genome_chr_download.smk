
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
        "logs/rule/analysis/{build}/log/download_chr{chr}.log",
    benchmark:
        "logs/rule/analysis/{build}/log/download_chr{chr}.benchmark.txt",
    wrapper:
        "0.69.0/bio/reference/ensembl-sequence"