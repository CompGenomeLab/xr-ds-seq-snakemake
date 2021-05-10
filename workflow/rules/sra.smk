rule sra_se:
    output:
        "resources/samples/{samples}.fastq.gz", 
    params:
        srr=lambda w: getSRR(w.samples, config["srr"]["codes"], 
            config["sample"]),
        name="{samples}",
    log:
        "logs/{samples}/{samples}_se_sra.log",
    benchmark:
        "logs/{samples}/{samples}_se_sra.benchmark.txt",
    wildcard_constraints:
        samples='|'.join([x for x in config["sample"]])
    conda:
        "../envs/sra.yaml"
    threads:
        8
    shell:
        """
        (echo "`date -R`: Downloading SRR files..." &&
        fasterq-dump \
        --threads {threads} \
        --progress {params.srr} \
        -t resources/samples/ \
        -o resources/samples/{params.name} &&
        gzip resources/samples/{params.name}*.fastq &&
        echo "`date -R`: Download is successful!" || 
        echo "`date -R`: Process failed...") \
        > {log} 2>&1
        """

rule sra_pe:
    output:
        "resources/samples/{samples}_1.fastq.gz", 
        "resources/samples/{samples}_2.fastq.gz", 
    params:
        srr=lambda w: getSRR(w.samples, config["srr"]["codes"], 
            config["sample"]),
        name="{samples}",
    log:
        "logs/{samples}/{samples}_pe_sra.log",
    benchmark:
        "logs/{samples}/{samples}_pe_sra.benchmark.txt",
    conda:
        "../envs/sra.yaml"
    threads:
        8
    shell:
        """
        (echo "`date -R`: Downloading SRR files..." &&
        fasterq-dump \
        --threads {threads} \
        --progress {params.srr} \
        -t resources/samples/ \
        -o resources/samples/{params.name} &&
        gzip resources/samples/{params.name}*.fastq &&
        echo "`date -R`: Download is successful!" || 
        echo "`date -R`: Process failed...") \
        > {log} 2>&1
        """
