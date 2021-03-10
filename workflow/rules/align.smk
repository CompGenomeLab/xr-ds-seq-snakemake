
rule bowtie2_se:
    input:
        sample=["results/{samples}/{samples}_cutadapt.fastq.gz"],
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2",
    output:
        "results/{samples}/{samples}_cutadapt_se_{build}.bam",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        extra=config["bowtie2"],
    threads: 4  
    log:
        "logs/{samples}/{samples}_{build}_bowtie2.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_bowtie2.benchmark.txt",
    conda:
        "../envs/align.yaml"
    shell:  
        """
        (bowtie2 --threads {threads} {params.extra} -x {params.ref_genome} \
        -U {input.sample[0]} | samtools view -Sbh -o {output} -) 2>{log}
        """

rule bowtie2_pe:
    input:
        sample=["results/{samples}/{samples}_cutadapt_1.fastq.gz", "results/{samples}/{samples}_cutadapt_2.fastq.gz"],
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2",
    output:
        "results/{samples}/{samples}_cutadapt_pe_{build}.bam",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        extra=config["bowtie2"],
    threads: 4  
    log:
        "logs/{samples}/{samples}_{build}_bowtie2.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_bowtie2.benchmark.txt",
    conda:
        "../envs/align.yaml"
    shell:  
        """
        (bowtie2 --threads {threads} {params.extra} -x {params.ref_genome} \
        -1 {input.sample[0]} -2 {input.sample[1]} | \
        samtools view -Sb -o {output} -) 2>{log}
        """