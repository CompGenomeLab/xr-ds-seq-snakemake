
rule genome_download:
    output:
        "resources/ref_genomes/{build}/genome_{build}.fa"
    params:
        link=config["genome"]["link"],
        build=config["genome"]["build"],
    log:
        "logs/rule/analysis/{build}/log/download.log",
    benchmark:
        "logs/rule/analysis/{build}/log/download.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Downloading {params.build} genome..." &&
        wget {params.link} -O {output}.gz &&
        gunzip {output}.gz &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """