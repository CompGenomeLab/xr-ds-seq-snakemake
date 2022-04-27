
rule genome_build:
    input:
        rules.genome_download.output,
    output:
        multiext(
        "resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
    log: 
        "logs/rule/analysis/{build}/log/bowtie2_build.log"
    benchmark:
        "logs/rule/analysis/{build}/log/bowtie2_build.benchmark.txt",
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
