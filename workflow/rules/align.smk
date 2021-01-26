
rule bowtie2_se:
    input:
        sample=["results/{dir}/{samples}{v}/{samples}_cutadapt.fastq.gz"],
    output:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_se.bam"
    params:
        ref_genome=getGenome,
        extra=lambda w: getSampleInfo(w, config["bowtie2"]),
    threads: 4  
    log:
        "results/{dir}/{samples}{v}/log/bowtie2.log"
    benchmark:
        "results/{dir}/{samples}{v}/log/bowtie2.benchmark.txt",
    conda:
        "../envs/align.yaml"
    shell:  
        """
        (bowtie2 --threads {threads} {params.extra} -x {params.ref_genome} \
        -U {input.sample[0]} | samtools view -Sbh -o {output} -) 2>{log}
        """

rule bowtie2_pe:
    input:
        sample=multiext(
            "results/{dir}/{samples}{v}/{samples}_cutadapt_",
            "1.fastq.gz", "2.fastq.gz"),
    output:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_pe.bam"
    params:
        ref_genome=getGenome,
        extra=lambda w: getSampleInfo(w, config["bowtie2"]),
    threads: 4  
    log:
        "results/{dir}/{samples}{v}/log/bowtie2.log"
    benchmark:
        "results/{dir}/{samples}{v}/log/bowtie2.benchmark.txt",
    conda:
        "../envs/align.yaml"
    shell:  
        """
        (bowtie2 --threads {threads} {params.extra} -x {params.ref_genome} \
        -1 {input.sample[0]} -2 {input.sample[1]} | \
        samtools view -Sbh -o {output} -) 2>{log}
        """