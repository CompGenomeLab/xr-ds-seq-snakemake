rule genome_variation:
    output:
        vcf="resources/ref_genomes/{build}/variation_{build}.vcf.gz",
    params:
        species=config["species"],
        build=config["build"],
        release=config["release"],
        type="all", # one of "all", "somatic", "structural_variation"
    log:
        "logs/rule/analysis/{build}/log/genome_variation.log",
    benchmark:
        "logs/rule/analysis/{build}/log/genome_variation.benchmark.txt",
    cache: True  
    wrapper:
        "0.69.0/bio/reference/ensembl-variation"
