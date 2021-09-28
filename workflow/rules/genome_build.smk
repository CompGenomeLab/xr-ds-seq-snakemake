
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
        {input.reference} \
        {params.name} &&
        echo "`date -R`: Success! Indexes are build." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """
