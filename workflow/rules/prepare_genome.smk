
rule genome_download:
    output:
        "resources/ref_genomes/{build}/genome_{build}.fa"
    params:
        link=config["genome"]["link"],
        build=config["genome"]["build"],
    log:
        "logs/rule/genome_download/{build}.log",
    benchmark:
        "logs/rule/genome_download/{build}.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Downloading {params.build} genome..." &&
        wget {params.link} -O {output}.gz &&
        gunzip {output}.gz &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """

rule genome_build:
    input:
        rules.genome_download.output,
    output:
        multiext(
        "resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    log: 
        "logs/rule/genome_build/{build}.log",
    benchmark:
        "logs/rule/genome_build/{build}.benchmark.txt",
    params:
        extra="", 
        name="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
    conda:
        "../envs/align.yaml" 
    threads: 
        16
    shell: 
        """
        (echo "`date -R`: Building indexes..." &&
        bowtie2-build --threads {threads} \
        {params.extra} \
        {input} \
        {params.name} &&
        echo "`date -R`: Success! Indexes are build." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """

rule genome_indexing:
    input:
        rules.genome_download.output,
    output:
        "resources/ref_genomes/{build}/genome_{build}.fa.fai",
    log: 
        "logs/rule/genome_indexing/{build}.log",
    benchmark:
        "logs/rule/genome_indexing/{build}.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        (echo "`date -R`: Creating fai file..." &&
        samtools faidx {input} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """

rule genome_index2ron:
    input:
        rules.genome_indexing.output,
    output:
        "resources/ref_genomes/{build}/genome_{build}.ron",
    log:
        "logs/rule/genome_index2ron/{build}.log",
    benchmark:
        "logs/rule/genome_index2ron/{build}.benchmark.txt",
    shell:
        "python3 workflow/scripts/idx2ron.py -i {input} -o {output} -l {log}"
