

rule bwa_index:
    input:
        rules.genome_download.output,
    output:
        multiext("resources/ref_genomes/{build}/genome_{build}.fa", 
        ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log: 
        "logs/prepare_genome/bwa_index/{build}.log",
    benchmark:
        "logs/prepare_genome/bwa_index/{build}.benchmark.txt",
    conda:
        "../envs/bwa.yaml"
    shell: 
        """    
        (echo "`date -R`: Building indexes..." &&
        bwa index {input} &&
        echo "`date -R`: Success! Indexes are build." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """

rule cutadapt_se:
    input:
        "resources/samples/{samples}.fastq.gz",
    output:
        fastq="results/{method}/{samples}/{samples}_cutadapt.fastq.gz",
        qc=report("results/{method}/{samples}/{samples}_cutadapt.qc.txt", category="QC"),  
    params:
        adapters=lambda w: getMethodParams(w, config["meta"], "adaptor", 
            config["XR"], config["DS"]),
        extra=lambda w: getMethodParams(w, config["meta"], "cutadapt", 
            config["XR"], config["DS"]),
    threads: 32
    log:
        "logs/rule/cutadapt_se/{samples}_{method}.log",
    benchmark:
        "logs/rule/cutadapt_se/{samples}_{method}.benchmark.txt",
    conda:
        "../envs/cutadapt.yaml"
    shell:
        """
        (echo "`date -R`: Trimming adapters..." &&
        cutadapt \
        -j {threads} \
        {params.adapters} \
        {params.extra} \
        -o {output.fastq} {input} \
        > {output.qc} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """

rule bowtie2_se:
    input:
        sample=[rules.cutadapt_se.output.fastq],
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2",
    output:
        sam=temp("results/{method}/{samples}/{samples}_cutadapt_se_{build}.sam"),
        bam="results/{method}/{samples}/{samples}_cutadapt_se_{build}.bam",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        extra="--seed 1 --reorder",
    threads: 32  
    log:
        report("logs/rule/bowtie2_se/{samples}_{build}_{method}.log", category="QC"),
    benchmark:
        "logs/rule/bowtie2_se/{samples}_{build}_{method}.benchmark.txt",
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

rule bowtie2_se_sensitive:
    input:
        sample=[rules.cutadapt_se.output.fastq],
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2",
    output:
        sam=temp("results/{method}/{samples}/{samples}_cutadapt_se_very_sensitive_{build}.sam"),
        bam="results/{method}/{samples}/{samples}_cutadapt_se_very_sensitive_{build}.bam",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        extra="--seed 1 --reorder --very-sensitive",
    threads: 32  
    log:
        report("logs/rule/bowtie2_se/{samples}_{build}_{method}_sensitive.log", category="QC"),
    benchmark:
        "logs/rule/bowtie2_se/{samples}_{build}_{method}_sensitive.benchmark.txt",
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

rule bwa_se:
    input:
        sample=[rules.cutadapt_se.output.fastq],
        genome=rules.genome_download.output,
        bwa=rules.bwa_index.output,
    output:
        sam=temp("results/{method}/{samples}/{samples}_cutadapt_se_bwa_{build}.sam"),
        sai=temp("results/{method}/{samples}/{samples}_cutadapt_se_bwa_{build}.sai"),
        bam="results/{method}/{samples}/{samples}_cutadapt_se_bwa_{build}.bam",
    params:
        "",
    threads: 32,
    log:
        "logs/rule/bwa_se/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/bwa_se/{samples}_{build}_{method}.benchmark.txt",
    conda:
        "../envs/bwa.yaml"
    shell:  
        """
        (echo "`date -R`: Aligning fastq files..." &&
        bwa aln -t {threads} {params} {input.genome} {input.sample} > {output.sai} &&
        bwa samse {input.genome} {output.sai} {input.sample} > {output.sam} &&
        echo "`date -R`: Success! Alignment is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1

        (echo "`date -R`: Converting sam to bam..." &&
        samtools sort -@ {threads} -o {output.bam} {output.sam} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; rm {output.bam}; exit 1; }}  )  >> {log} 2>&1
        """

rule bwa_se_sensitive1:
    input:
        sample=[rules.cutadapt_se.output.fastq],
        genome=rules.genome_download.output,
        bwa=rules.bwa_index.output,
    output:
        sam=temp("results/{method}/{samples}/{samples}_cutadapt_se_bwa_sensitive1_{build}.sam"),
        sai=temp("results/{method}/{samples}/{samples}_cutadapt_se_bwa_sensitive1_{build}.sai"),
        bam="results/{method}/{samples}/{samples}_cutadapt_se_bwa_sensitive1_{build}.bam",
    params:
        "-l 15",
    threads: 32,
    log:
        "logs/rule/bwa_se/{samples}_{build}_{method}_sensitive1.log",
    benchmark:
        "logs/rule/bwa_se/{samples}_{build}_{method}_sensitive1.benchmark.txt",
    conda:
        "../envs/bwa.yaml"
    shell:  
        """
        (echo "`date -R`: Aligning fastq files..." &&
        bwa aln -t {threads} {params} {input.genome} {input.sample} > {output.sai} &&
        bwa samse {input.genome} {output.sai} {input.sample} > {output.sam} &&
        echo "`date -R`: Success! Alignment is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1

        (echo "`date -R`: Converting sam to bam..." &&
        samtools sort -@ {threads} -o {output.bam} {output.sam} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; rm {output.bam}; exit 1; }}  )  >> {log} 2>&1
        """

rule bwa_se_sensitive2:
    input:
        sample=[rules.cutadapt_se.output.fastq],
        genome=rules.genome_download.output,
        bwa=rules.bwa_index.output,
    output:
        sam=temp("results/{method}/{samples}/{samples}_cutadapt_se_bwa_sensitive2_{build}.sam"),
        sai=temp("results/{method}/{samples}/{samples}_cutadapt_se_bwa_sensitive2_{build}.sai"),
        bam="results/{method}/{samples}/{samples}_cutadapt_se_bwa_sensitive2_{build}.bam",
    params:
        "-n 3",
    threads: 32,
    log:
        "logs/rule/bwa_se/{samples}_{build}_{method}_sensitive2.log",
    benchmark:
        "logs/rule/bwa_se/{samples}_{build}_{method}_sensitive2.benchmark.txt",
    conda:
        "../envs/bwa.yaml"
    shell:  
        """
        (echo "`date -R`: Aligning fastq files..." &&
        bwa aln -t {threads} {params} {input.genome} {input.sample} > {output.sai} &&
        bwa samse {input.genome} {output.sai} {input.sample} > {output.sam} &&
        echo "`date -R`: Success! Alignment is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1

        (echo "`date -R`: Converting sam to bam..." &&
        samtools sort -@ {threads} -o {output.bam} {output.sam} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; rm {output.bam}; exit 1; }}  )  >> {log} 2>&1
        """

rule bwa_se_sensitive3:
    input:
        sample=[rules.cutadapt_se.output.fastq],
        genome=rules.genome_download.output,
        bwa=rules.bwa_index.output,
    output:
        sam=temp("results/{method}/{samples}/{samples}_cutadapt_se_bwa_sensitive3_{build}.sam"),
        sai=temp("results/{method}/{samples}/{samples}_cutadapt_se_bwa_sensitive3_{build}.sai"),
        bam="results/{method}/{samples}/{samples}_cutadapt_se_bwa_sensitive3_{build}.bam",
    params:
        "-e 0",
    threads: 32,
    log:
        "logs/rule/bwa_se/{samples}_{build}_{method}_sensitive3.log",
    benchmark:
        "logs/rule/bwa_se/{samples}_{build}_{method}_sensitive3.benchmark.txt",
    conda:
        "../envs/bwa.yaml"
    shell:  
        """
        (echo "`date -R`: Aligning fastq files..." &&
        bwa aln -t {threads} {params} {input.genome} {input.sample} > {output.sai} &&
        bwa samse {input.genome} {output.sai} {input.sample} > {output.sam} &&
        echo "`date -R`: Success! Alignment is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1

        (echo "`date -R`: Converting sam to bam..." &&
        samtools sort -@ {threads} -o {output.bam} {output.sam} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; rm {output.bam}; exit 1; }}  )  >> {log} 2>&1
        """

rule bwa_se_sensitive4:
    input:
        sample=[rules.cutadapt_se.output.fastq],
        genome=rules.genome_download.output,
        bwa=rules.bwa_index.output,
    output:
        sam=temp("results/{method}/{samples}/{samples}_cutadapt_se_bwa_sensitive4_{build}.sam"),
        sai=temp("results/{method}/{samples}/{samples}_cutadapt_se_bwa_sensitive4_{build}.sai"),
        bam="results/{method}/{samples}/{samples}_cutadapt_se_bwa_sensitive4_{build}.bam",
    params:
        "-k 3",
    threads: 32,
    log:
        "logs/rule/bwa_se/{samples}_{build}_{method}_sensitive4.log",
    benchmark:
        "logs/rule/bwa_se/{samples}_{build}_{method}_sensitive4.benchmark.txt",
    conda:
        "../envs/bwa.yaml"
    shell:  
        """
        (echo "`date -R`: Aligning fastq files..." &&
        bwa aln -t {threads} {params} {input.genome} {input.sample} > {output.sai} &&
        bwa samse {input.genome} {output.sai} {input.sample} > {output.sam} &&
        echo "`date -R`: Success! Alignment is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1

        (echo "`date -R`: Converting sam to bam..." &&
        samtools sort -@ {threads} -o {output.bam} {output.sam} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; rm {output.bam}; exit 1; }}  )  >> {log} 2>&1
        """

rule bwa_se_sensitive5:
    input:
        sample=[rules.cutadapt_se.output.fastq],
        genome=rules.genome_download.output,
        bwa=rules.bwa_index.output,
    output:
        sam=temp("results/{method}/{samples}/{samples}_cutadapt_se_bwa_sensitive5_{build}.sam"),
        sai=temp("results/{method}/{samples}/{samples}_cutadapt_se_bwa_sensitive5_{build}.sai"),
        bam="results/{method}/{samples}/{samples}_cutadapt_se_bwa_sensitive5_{build}.bam",
    params:
        "-l 15 -n 3 -e 0 -k 3",
    threads: 32,
    log:
        "logs/rule/bwa_se/{samples}_{build}_{method}_sensitive5.log",
    benchmark:
        "logs/rule/bwa_se/{samples}_{build}_{method}_sensitive5.benchmark.txt",
    conda:
        "../envs/bwa.yaml"
    shell:  
        """
        (echo "`date -R`: Aligning fastq files..." &&
        bwa aln -t {threads} {params} {input.genome} {input.sample} > {output.sai} &&
        bwa samse {input.genome} {output.sai} {input.sample} > {output.sam} &&
        echo "`date -R`: Success! Alignment is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1

        (echo "`date -R`: Converting sam to bam..." &&
        samtools sort -@ {threads} -o {output.bam} {output.sam} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; rm {output.bam}; exit 1; }}  )  >> {log} 2>&1
        """