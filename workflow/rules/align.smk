
rule bowtie2_se:
    input:
        sample=["results/{samples}/{samples}_cutadapt.fastq.gz"],
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2",
    output:
        sam=temp("results/{samples}/{samples}_cutadapt_se_{build}.sam"),
        bam="results/{samples}/{samples}_cutadapt_se_{build}.bam",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        extra="",
    threads: 16  
    log:
        "logs/{samples}/{samples}_{build}_bowtie2.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_bowtie2.benchmark.txt",
    conda:
        "../envs/align.yaml"
    shell:  
        """
        (echo "`date -R`: Aligning fastq file..." &&
        bowtie2 \
        --threads {threads} \
        {params.extra} \
        -x {params.ref_genome} \
        -U {input.sample[0]} -S {output.sam} &&
        echo "`date -R`: Success! Alignment is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1

        (echo "`date -R`: Converting sam to bam..." &&
        samtools view -Sbh -o {output.bam} {output.sam} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; rm {output.bam}; exit 1; }}  )  >> {log} 2>&1
        """

rule bowtie2_pe:
    input:
        sample=["results/{samples}/{samples}_cutadapt_1.fastq.gz", "results/{samples}/{samples}_cutadapt_2.fastq.gz"],
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2",
    output:
        sam=temp("results/{samples}/{samples}_cutadapt_pe_{build}.sam"),
        bam="results/{samples}/{samples}_cutadapt_pe_{build}.bam",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        extra="-X 1000",
    threads: 16  
    log:
        "logs/{samples}/{samples}_{build}_bowtie2.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_bowtie2.benchmark.txt",
    conda:
        "../envs/align.yaml"
    shell:  
        """
        (echo "`date -R`: Aligning fastq files..." &&
        bowtie2 \
        --threads {threads} \
        {params.extra} \
        -x {params.ref_genome} \
        -1 {input.sample[0]} -2 {input.sample[1]} -S {output.sam} &&
        echo "`date -R`: Success! Alignment is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1

        (echo "`date -R`: Converting sam to bam..." &&
        samtools view -Sb -o {output.bam} {output.sam} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; rm {output.bam}; exit 1; }}  )  >> {log} 2>&1
        """