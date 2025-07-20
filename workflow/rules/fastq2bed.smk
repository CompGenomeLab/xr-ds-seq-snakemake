rule cutadapt_se:
    input:
        "resources/samples/{samples}.fastq.gz",
    output:
        fastq=temp("results/{method}/{samples}/{samples}_cutadapt.fastq.gz"),
        qc=report("results/{method}/{samples}/{samples}_cutadapt.qc.txt", category="QC"),  
    params:
        adapters=lambda w: getMethodParams(w, config["meta"], "adaptor", 
            config["XR"], config["DS"]),
        extra=lambda w: getMethodParams(w, config["meta"], "cutadapt", 
            config["XR"], config["DS"]),
    threads: 16
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

rule cutadapt_pe:
    input:
        fq1="resources/samples/{samples}_1.fastq.gz", 
        fq2="resources/samples/{samples}_2.fastq.gz", 
    output:
        fastq1=temp("results/{method}/{samples}/{samples}_cutadapt_1.fastq.gz"),
        fastq2=temp("results/{method}/{samples}/{samples}_cutadapt_2.fastq.gz"),
        qc=report("results/{method}/{samples}/{samples}_cutadapt.qc.txt", category="QC"),
    params:
        adapters=lambda w: getMethodParams(w, config["meta"], "adaptor", 
            config["XR"], config["DS"]),
        extra=lambda w: getMethodParams(w, config["meta"], "cutadapt", 
            config["XR"], config["DS"]),
    threads: 16
    log:
        "logs/rule/cutadapt_pe/{samples}_{method}.log",
    benchmark:
        "logs/rule/cutadapt_pe/{samples}_{method}.benchmark.txt",
    conda:
        "../envs/cutadapt.yaml"
    shell:
        """
        (echo "`date -R`: Trimming adapters for paired-end..." &&
        cutadapt \
        -j {threads} \
        {params.adapters} \
        {params.extra} \
        -o {output.fastq1} \
        -p {output.fastq2} \
        {input} \
        > {output.qc} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """

rule process_fastq_T2C:
    input:
        "results/{method}/{samples}/{samples}_cutadapt_se_{build}_unmapped_filtered_overrep.fastq"
    output:
        "results/{method}/{samples}/{samples}_T2C_{build}_cutadapt.fastq.gz"
    log:
        "logs/rule/process_fastq_T2C/{samples}_{build}_{method}.log",
    shell:
        "bash workflow/scripts/t2c_damsite.sh {input} {output} > {log} 2>&1"

rule mask_positions:
    input:
        "results/{method}/{samples}/{samples}_cutadapt_se_{build}_unmapped_filtered_overrep.fastq"
    output:
        "results/{method}/{samples}/{samples}_masked_{build}_cutadapt.fastq.gz"
    log:
        "logs/rule/mask_positions/{samples}_{build}_{method}.log",
    shell:
        """
        ( paste - - - - < {input} | \
        awk -F'\\t' '{{if($2!="") {{ \
            len=length($2); \
            $2=substr($2,1,len-8) "NN" substr($2,len-5); \
            if(length($2)==len) {{ \
                print $1"\\n"$2"\\n"$3"\\n"$4 \
            }} \
        }}}}' | \
        gzip > {output} ) > {log} 2>&1
        """

rule bowtie2_se_filtered:
    input:
        sample=["results/{method}/{samples}/{samples}_cutadapt_se_filtered_overrep.fastq"],
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2",
    output:
        sam=temp("results/{method}/{samples}/{samples}_filt_cutadapt_se_{build}.sam"),
        bam="results/{method}/{samples}/{samples}_filt_cutadapt_se_{build}.bam",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        extra="--seed 1 --reorder",
    threads: 16  
    log:
        report("logs/rule/bowtie2_se_filtered/{samples}_{build}_{method}.log", category="QC"),
    benchmark:
        "logs/rule/bowtie2_se_filtered/{samples}_{build}_{method}.benchmark.txt",
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

rule bowtie2_se_T2C:
    input:
        sample=[rules.process_fastq_T2C.output],
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2",
    output:
        sam=temp("results/{method}/{samples}/{samples}_T2C_cutadapt_se_{build}.sam"),
        bam="results/{method}/{samples}/{samples}_T2C_cutadapt_se_{build}.bam",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        extra="--seed 1 --reorder",
    threads: 16  
    log:
        report("logs/rule/bowtie2_se_T2C/{samples}_{build}_{method}.log", category="QC"),
    benchmark:
        "logs/rule/bowtie2_se_T2C/{samples}_{build}_{method}.benchmark.txt",
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

rule bowtie2_se_masked:
    input:
        sample=[rules.mask_positions.output],
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2",
    output:
        sam=temp("results/{method}/{samples}/{samples}_masked_cutadapt_se_{build}.sam"),
        bam="results/{method}/{samples}/{samples}_masked_cutadapt_se_{build}.bam",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        extra="--seed 1 --reorder",
    threads: 16  
    log:
        report("logs/rule/bowtie2_se_masked/{samples}_{build}_{method}.log", category="QC"),
    benchmark:
        "logs/rule/bowtie2_se_masked/{samples}_{build}_{method}.benchmark.txt",
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

rule bowtie2_se:
    input:
        sample=[rules.cutadapt_se.output.fastq],
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2",
    output:
        sam=temp("results/{method}/{samples}/{samples}_cutadapt_se_{build}.sam"),
        bam="results/{method}/{samples}/{samples}_cutadapt_se_{build}.bam",
        unmapped="results/{method}/{samples}/{samples}_cutadapt_se_{build}_unmapped.fastq",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        extra="--seed 1 --reorder",
    threads: 16  
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
        --un {output.unmapped} \
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
        sample=[rules.cutadapt_pe.output.fastq1, rules.cutadapt_pe.output.fastq2],
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2",
    output:
        sam=temp("results/{method}/{samples}/{samples}_cutadapt_pe_{build}.sam"),
        bam="results/{method}/{samples}/{samples}_cutadapt_pe_{build}.bam",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        extra="-X 1000 --seed 1 --reorder",
    threads: 16  
    log:
        report("logs/rule/bowtie2_pe/{samples}_{build}_{method}.log", category="QC"),
    benchmark:
        "logs/rule/bowtie2_pe/{samples}_{build}_{method}.benchmark.txt",
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

rule bowtie2_se_input:
    input:
        sample=["resources/input/{samples}.fastq.gz"],
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2",
    output:
        sam=temp("results/input/{samples}/{samples}_se_{build}.sam"),
        bam="results/input/{samples}/{samples}_se_{build}.bam",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        extra="--seed 1 --reorder",
    threads: 16  
    log:
        report("logs/rule/bowtie2_se_input/{samples}_{build}.log", category="QC"),
    benchmark:
        "logs/rule/bowtie2_se_input/{samples}_{build}.benchmark.txt",
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

rule bowtie2_pe_input:
    input:
        sample=["resources/input/{samples}_1.fastq.gz", "resources/input/{samples}_2.fastq.gz"],
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2",
    output:
        sam=temp("results/input/{samples}/{samples}_pe_{build}.sam"),
        bam="results/input/{samples}/{samples}_pe_{build}.bam",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        extra="-X 1000 --seed 1 --reorder",
    threads: 16  
    log:
        report("logs/rule/bowtie2_pe_input/{samples}_{build}.log", category="QC"),
    benchmark:
        "logs/rule/bowtie2_pe_input/{samples}_{build}.benchmark.txt",
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

rule sort_rmPG_se_filtered:
    input:
        rules.bowtie2_se_filtered.output.bam,
    output:
        header=temp("results/{method}/{samples}/{samples}_cutadapt_se_{build}_filtered_header.txt"),
        sort=temp("results/{method}/{samples}/{samples}_cutadapt_se_{build}_filtered_samSorted.bam"),
    params:
        tmpdir="results/{method}/{samples}/",
    threads: 16
    log:
        "logs/rule/sort_rmPG_se_filtered/{samples}_{method}_{build}.log",
    benchmark:
        "logs/rule/sort_rmPG_se_filtered/{samples}_{method}_{build}.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        (echo "`date -R`: Parsing the headers..." &&
        samtools view -H {input} | grep -v "^@PG" > {output.header} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.header}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Parsing the headers..." &&
        samtools reheader {output.header} {input} | 
        samtools sort -o {output.sort} -@ {threads} -T {params.tmpdir} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.sort}; exit 1; }}  ) >> {log} 2>&1
        """

rule sort_rmPG_se:
    input:
        rules.bowtie2_se.output.bam,
    output:
        header=temp("results/{method}/{samples}/{samples}_cutadapt_se_{build}_header.txt"),
        sort=temp("results/{method}/{samples}/{samples}_cutadapt_se_{build}_samSorted.bam"),
    params:
        tmpdir="results/{method}/{samples}/",
    threads: 16
    log:
        "logs/rule/sort_rmPG_se/{samples}_{method}_{build}.log",
    benchmark:
        "logs/rule/sort_rmPG_se/{samples}_{method}_{build}.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        (echo "`date -R`: Parsing the headers..." &&
        samtools view -H {input} | grep -v "^@PG" > {output.header} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.header}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Parsing the headers..." &&
        samtools reheader {output.header} {input} | 
        samtools sort -o {output.sort} -@ {threads} -T {params.tmpdir} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.sort}; exit 1; }}  ) >> {log} 2>&1
        """

rule sort_rmPG_pe:
    input:
        rules.bowtie2_pe.output.bam,
    output:
        header=temp("results/{method}/{samples}/{samples}_cutadapt_pe_{build}_header.txt"),
        sort=temp("results/{method}/{samples}/{samples}_cutadapt_pe_{build}_samSorted.bam"),
    params:
        tmpdir="results/{method}/{samples}/",
    threads: 16
    log:
        "logs/rule/sort_rmPG_pe/{samples}_{method}_{build}.log",
    benchmark:
        "logs/rule/sort_rmPG_pe/{samples}_{method}_{build}.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        (echo "`date -R`: Parsing the headers excluding the @PG headers..." &&
        samtools view -H {input} | grep -v "^@PG" > {output.header} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.header}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Reheader and sort..." &&
        samtools reheader {output.header} {input} | 
        samtools sort -o {output.sort} -@ {threads} -T {params.tmpdir} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.sort}; exit 1; }}  ) >> {log} 2>&1
        """

rule mark_duplicates_se_filtered:
    input:
        rules.sort_rmPG_se_filtered.output.sort,
    output:
        bam="results/{method}/{samples}/{samples}_dedup_cutadapt_se_filtered_{build}.bam",
        metrics="results/{method}/{samples}/{samples}_cutadapt_se_filtered_dedup_{build}.metrics.txt",
    log:
        "logs/rule/mark_duplicates_se_filtered/{samples}_{method}_{build}.log",
    params:
        extra="--REMOVE_DUPLICATES true",
        tmpdir="results/{method}/{samples}/",
    conda:
        "../envs/picard.yaml"
    shell:
        """
        (echo "`date -R`: Removing duplicates..." &&
        picard MarkDuplicates \
        {params.extra} \
        --INPUT {input} \
        --TMP_DIR {params.tmpdir} \
        --OUTPUT {output.bam} \
        --METRICS_FILE {output.metrics} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """

rule mark_duplicates_se:
    input:
        rules.sort_rmPG_se.output.sort,
    output:
        bam="results/{method}/{samples}/{samples}_dedup_cutadapt_se_{build}.bam",
        metrics="results/{method}/{samples}/{samples}_cutadapt_se_dedup_{build}.metrics.txt",
    log:
        "logs/rule/mark_duplicates_se/{samples}_{method}_{build}.log",
    params:
        extra="--REMOVE_DUPLICATES true",
        tmpdir="results/{method}/{samples}/",
    conda:
        "../envs/picard.yaml"
    shell:
        """
        (echo "`date -R`: Removing duplicates..." &&
        picard MarkDuplicates \
        {params.extra} \
        --INPUT {input} \
        --TMP_DIR {params.tmpdir} \
        --OUTPUT {output.bam} \
        --METRICS_FILE {output.metrics} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """

rule mark_duplicates_pe:
    input:
        rules.sort_rmPG_pe.output.sort,
    output:
        bam="results/{method}/{samples}/{samples}_dedup_cutadapt_pe_{build}.bam",
        metrics="results/{method}/{samples}/{samples}_cutadapt_pe_dedup_{build}.metrics.txt",
    log:
        "logs/rule/mark_duplicates_pe/{samples}_{method}_{build}.log",
    params:
        extra="--REMOVE_DUPLICATES true",
        tmpdir="results/{method}/{samples}/",
    conda:
        "../envs/picard.yaml"
    shell:
        """
        (echo "`date -R`: Removing duplicates..." &&
        picard MarkDuplicates \
        {params.extra} \
        --INPUT {input} \
        --TMP_DIR {params.tmpdir} \
        --OUTPUT {output.bam} \
        --METRICS_FILE {output.metrics} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """

rule bam2bed_se_filtered:
    input:
        rules.mark_duplicates_se_filtered.output.bam,
    output:
        bed="results/{method}/{samples}/{samples}_{build}_se_filtered.bed",
        bam="results/{method}/{samples}/{samples}_{build}_se_filtered_sortedbyCoordinates.bam",
        idx="results/{method}/{samples}/{samples}_{build}_se_filtered_sortedbyCoordinates.bam.bai",
    params:
        q_trim=lambda w: getMethodParams(w, config["meta"], "samtools", 
            config["XR"], config["DS"]), 
        tmpdir="results/{method}/{samples}/",
    threads: 16
    log:
        "logs/rule/bam2bed_se_filtered/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/bam2bed_se_filtered/{samples}_{build}_{method}.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Sorting (coordinates) bam file..." &&
        samtools sort {input} -o {output.bam} -@ {threads} -T {params.tmpdir} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Index bam file..." &&
        samtools index {output.bam} {output.idx} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Processing bam file..." && 
        samtools view {params.q_trim} -b {input} |&
        bedtools bamtobed > {output.bed} &&
        echo "`date -R`: Success! Bam file converted to bed format." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1
        """

rule bam2bed_se:
    input:
        rules.mark_duplicates_se.output.bam,
    output:
        bed="results/{method}/{samples}/{samples}_{build}_se.bed",
        bam="results/{method}/{samples}/{samples}_{build}_se_sortedbyCoordinates.bam",
        idx="results/{method}/{samples}/{samples}_{build}_se_sortedbyCoordinates.bam.bai",
    params:
        q_trim=lambda w: getMethodParams(w, config["meta"], "samtools", 
            config["XR"], config["DS"]), 
        tmpdir="results/{method}/{samples}/",
    threads: 16
    log:
        "logs/rule/bam2bed_se/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/bam2bed_se/{samples}_{build}_{method}.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Sorting (coordinates) bam file..." &&
        samtools sort {input} -o {output.bam} -@ {threads} -T {params.tmpdir} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Index bam file..." &&
        samtools index {output.bam} {output.idx} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Processing bam file..." && 
        samtools view {params.q_trim} -b {input} |&
        bedtools bamtobed > {output.bed} &&
        echo "`date -R`: Success! Bam file converted to bed format." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1
        """

rule bam2bed_pe:
    input:
        rules.mark_duplicates_pe.output.bam,
    output:
        bed="results/{method}/{samples}/{samples}_{build}_pe.bed",
        bam=temp("results/{method}/{samples}/{samples}_{build}_sorted.bam"),
        bam2="results/{method}/{samples}/{samples}_{build}_pe_sortedbyCoordinates.bam",
        idx="results/{method}/{samples}/{samples}_{build}_pe_sortedbyCoordinates.bam.bai",
    params:
        q_trim=lambda w: getMethodParams(w, config["meta"], "samtools", 
            config["XR"], config["DS"]), 
        tmpdir="results/{method}/{samples}/",
    threads: 16
    log:
        "logs/rule/bam2bed_pe/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/bam2bed_pe/{samples}_{build}_{method}.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Sorting (name) bam file..." &&
        samtools sort -n {input} -o {output.bam} \
        -@ {threads} -T {params.tmpdir} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Sorting (coordinates) bam file..." &&
        samtools sort {input} -o {output.bam2} \
        -@ {threads} -T {params.tmpdir} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Index bam file..." &&
        samtools index {output.bam2} {output.idx} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Processing bam file..." &&
        samtools view {params.q_trim} {output.bam} |&
        bedtools bamtobed -bedpe -mate1 |&
        awk '{{\
            if ($9=="+")\
                print $1"\\t"$2"\\t"$6"\\t"$7"\\t"$8"\\t"$9;\
            else if ($9=="-")\
                print $1"\\t"$5"\\t"$3"\\t"$7"\\t"$8"\\t"$9;\
            }}' > {output.bed} &&
        echo "`date -R`: Success! Bam file converted to bed format." || 
        {{ echo "`date -R`: Process failed..."; rm {output.bed}; exit 1; }}  ) >> {log} 2>&1
        """

rule bam2bed_se_input:
    input:
        rules.bowtie2_se_input.output.bam,
    output:
        bed="results/input/{samples}/{samples}_{build}_se.bed",
        bam="results/input/{samples}/{samples}_{build}_se_sortedbyCoordinates.bam",
        idx="results/input/{samples}/{samples}_{build}_se_sortedbyCoordinates.bam.bai",
    params:
        q_trim="-q 20",
        tmpdir="results/input/{samples}/",
    threads: 16
    log:
        "logs/rule/bam2bed_se_input/{samples}_{build}.log",
    benchmark:
        "logs/rule/bam2bed_se_input/{samples}_{build}.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Sorting (coordinates) bam file..." &&
        samtools sort {input} -o {output.bam} -@ {threads} -T {params.tmpdir} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Index bam file..." &&
        samtools index {output.bam} {output.idx} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Processing bam file..." && 
        samtools view {params.q_trim} -b {input} |&
        bedtools bamtobed > {output.bed} &&
        echo "`date -R`: Success! Bam file converted to bed format." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1
        """

rule bam2bed_pe_input:
    input:
        rules.bowtie2_pe_input.output.bam,
    output:
        bed="results/input/{samples}/{samples}_{build}_pe.bed",
        bam=temp("results/input/{samples}/{samples}_{build}_sorted.bam"),
        bam2="results/input/{samples}/{samples}_{build}_pe_sortedbyCoordinates.bam",
        idx="results/input/{samples}/{samples}_{build}_pe_sortedbyCoordinates.bam.bai",
    params:
        q_trim="-q 20 -bf 0x2",
        tmpdir="results/input/{samples}/",
    threads: 16
    log:
        "logs/rule/bam2bed_pe_input/{samples}_{build}.log",
    benchmark:
        "logs/rule/bam2bed_pe_input/{samples}_{build}.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Sorting (name) bam file..." &&
        samtools sort -n {input} -o {output.bam} \
        -@ {threads} -T {params.tmpdir} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Sorting (coordinates) bam file..." &&
        samtools sort {input} -o {output.bam2} \
        -@ {threads} -T {params.tmpdir} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Index bam file..." &&
        samtools index {output.bam2} {output.idx} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Processing bam file..." &&
        samtools view {params.q_trim} {output.bam} |&
        bedtools bamtobed -bedpe -mate1 |&
        awk '{{\
            if ($9=="+")\
                print $1"\\t"$2"\\t"$6"\\t"$7"\\t"$8"\\t"$9;\
            else if ($9=="-")\
                print $1"\\t"$5"\\t"$3"\\t"$7"\\t"$8"\\t"$9;\
            }}' > {output.bed} &&
        echo "`date -R`: Success! Bam file converted to bed format." || 
        {{ echo "`date -R`: Process failed..."; rm {output.bed}; exit 1; }}  ) >> {log} 2>&1
        """

rule filter_adapter_reads:
    input:
        "results/{method}/{samples}/{samples}_cutadapt_se_{build}_unmapped.fastq"
    output:
        "results/{method}/{samples}/{samples}_cutadapt_se_{build}_unmapped_filtered.fastq"
    log:
        "logs/rule/filter_adapter_reads/{samples}_{build}_{method}.log"
    shell:
        """
        (echo "`date -R`: Filtering out adapter reads..." &&
        paste - - - - < {input} | \
        awk -F'\\t' '{{if(index($2, "TGGAATTCTCGGGTGCCAA") == 0) {{print $1"\\n"$2"\\n"$3"\\n"$4}}}}' > {output} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """

rule extract_overrepresented_unmapped:
    input:
        zip="results/{method}/{samples}/{samples}_cutadapt_se_{build}_unmapped_fastqc.zip"
    output:
        sequences="results/{method}/{samples}/{samples}_cutadapt_se_{build}_unmapped_overrepresented.txt"
    log:
        "logs/rule/extract_overrepresented_unmapped/{samples}_{build}_{method}.log"
    shell:
        """
        unzip -p {input.zip} */fastqc_data.txt | \
        awk '/^>>Overrepresented sequences/,/^>>END_MODULE/ {{if($0 !~ /^>>/) print $1}}' > {output.sequences}
        """

rule extract_overrepresented:
    input:
        zip="results/processed_files/{samples}_cutadapt_fastqc.zip"
    output:
        sequences="results/{method}/{samples}/{samples}_cutadapt_se_overrepresented.txt"
    log:
        "logs/rule/extract_overrepresented/{samples}_{method}.log"
    shell:
        """
        unzip -p {input.zip} */fastqc_data.txt | \
        awk '/^>>Overrepresented sequences/,/^>>END_MODULE/ {{if($0 !~ /^>>/) print $1}}' > {output.sequences}
        """

rule filter_overrepresented:
    input:
        fastq="results/{method}/{samples}/{samples}_cutadapt.fastq.gz",
        sequences="results/{method}/{samples}/{samples}_cutadapt_se_overrepresented.txt"
    output:
        "results/{method}/{samples}/{samples}_cutadapt_se_filtered_overrep.fastq"
    log:
        "logs/rule/filter_overrepresented/{samples}_{method}.log"
    shell:
        """
        (echo "`date -R`: Filtering out reads that exactly match overrepresented sequences..." &&
        gunzip -c {input.fastq} | \
        paste - - - - | \
        awk -F'\\t' -v seqfile={input.sequences} '
        BEGIN {{while((getline seq < seqfile) > 0) seqs[seq]=1}}
        {{keep=1; \
            for(s in seqs) {{ \
                if($2 == s) {{ \
                    keep=0; \
                    break \
                }} \
            }}; \
            if(keep) print $1"\\n"$2"\\n"$3"\\n"$4 \
        }}' > {output} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """

rule filter_overrepresented_unmapped:
    input:
        fastq="results/{method}/{samples}/{samples}_cutadapt_se_{build}_unmapped.fastq",
        sequences="results/{method}/{samples}/{samples}_cutadapt_se_{build}_unmapped_overrepresented.txt"
    output:
        "results/{method}/{samples}/{samples}_cutadapt_se_{build}_unmapped_filtered_overrep.fastq"
    log:
        "logs/rule/filter_overrepresented_unmapped/{samples}_{build}_{method}.log"
    shell:
        """
        (echo "`date -R`: Filtering out reads that exactly match overrepresented sequences..." &&
        paste - - - - < {input.fastq} | \
        awk -F'\\t' -v seqfile={input.sequences} '
        BEGIN {{while((getline seq < seqfile) > 0) seqs[seq]=1}}
        {{keep=1; \
            for(s in seqs) {{ \
                if($2 == s) {{ \
                    keep=0; \
                    break \
                }} \
            }}; \
            if(keep) print $1"\\n"$2"\\n"$3"\\n"$4 \
        }}' > {output} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """
