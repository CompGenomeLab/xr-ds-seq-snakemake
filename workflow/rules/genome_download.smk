
rule get_genome:
    output:
        "resources/ref_genomes/{build}/genome_{build}.fa"
    params:
        species=lambda w: getGenomeInfo(w, config["species"]),
        datatype=lambda w: getGenomeInfo(w, config["datatype"]),
        build=lambda w: getGenomeInfo(w, config["build"]),
        release=lambda w: getGenomeInfo(w, config["release"]),
    log:
        "resources/ref_genomes/{build}/log/download.log",
    benchmark:
        "resources/ref_genomes/{build}/log/download.benchmark.txt",
    wrapper:
        "0.69.0/bio/reference/ensembl-sequence"
