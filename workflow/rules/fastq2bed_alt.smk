
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
        {{
            keep=1;
            for(s in seqs) {{
                if($2 == s) {{
                    keep=0;
                    break
                }}
            }}
            if(keep && length($2) > 0)
                print $1"\\n"$2"\\n"$3"\\n"$4
        }}' > {output} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """

rule fastx_collapser_se_filtered:
    input:
        fastq=rules.filter_overrepresented.output,
    output:
        fasta="results/{method}/{samples}/{samples}_cutadapt_se_filtered_overrep_collapsed.fasta.gz",
    log:
        "logs/rule/fastx_collapser_se_filtered/{samples}_{method}.log",
    benchmark:
        "logs/rule/fastx_collapser_se_filtered/{samples}_{method}.benchmark.txt",
    conda:
        "../envs/fastx_collapser.yaml"
    shell:
        """
        (echo "`date -R`: Collapsing duplicate reads (filtered FASTQ)..." &&
        cat {input.fastq} | \
        fastx_collapser -v | \
        gzip > {output.fasta} &&
        echo "`date -R`: Success! Reads are collapsed." || 
        {{ echo "`date -R`: Process failed..."; rm -f {output.fasta}; exit 1; }}  ) > {log} 2>&1
        """

rule bowtie2_se_filtered:
    input:
        sample=["results/{method}/{samples}/{samples}_cutadapt_se_filtered_overrep_collapsed.fasta.gz"],
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2",
    output:
        sam=temp("results/{method}/{samples}/{samples}_filt_cutadapt_se_{build}.sam"),
        bam="results/{method}/{samples}/{samples}_filt_cutadapt_se_{build}.bam",
        unmapped="results/{method}/{samples}/{samples}_cutadapt_se_{build}_filtered_overrep_collapsed_unmapped.fasta",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        extra="--seed 1 --reorder -f",
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
        --un {output.unmapped} \
        -U {input.sample[0]} -S {output.sam} &&
        echo "`date -R`: Success! Alignment is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1

        (echo "`date -R`: Converting sam to bam..." &&
        samtools view -Sbh -o {output.bam} {output.sam} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; rm {output.bam}; exit 1; }}  )  >> {log} 2>&1
        """

rule process_fasta_T2C:
    input:
        "results/{method}/{samples}/{samples}_cutadapt_se_{build}_filtered_overrep_collapsed_unmapped.fasta"
    output:
        "results/{method}/{samples}/{samples}_T2C_{build}_cutadapt.fasta.gz"
    log:
        "logs/rule/process_fastq_T2C/{samples}_{build}_{method}.log",
    shell:
        """
        (echo "`date -R`: Converting TT to CC across reads..." &&
        python workflow/scripts/tt_to_cc_converter.py {input} {output} &&
        echo "`date -R`: Success!") > {log} 2>&1
        """

rule bowtie2_se_T2C:
    input:
        sample=[rules.process_fasta_T2C.output],
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2",
    output:
        sam=temp("results/{method}/{samples}/{samples}_T2C_cutadapt_se_{build}.sam"),
        bam="results/{method}/{samples}/{samples}_T2C_cutadapt_se_{build}.bam",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        extra="--seed 1 --reorder -f",
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

rule merge_bams:
    input:
        bam1=rules.bowtie2_se_filtered.output.bam,
        bam2=rules.bowtie2_se_T2C.output.bam,
    output:
        merged_bam="results/{method}/{samples}/{samples}_{build}_merged.bam",
        merged_bam_idx="results/{method}/{samples}/{samples}_{build}_merged.bam.bai",
    params:
        tmpdir="results/{method}/{samples}/",
    threads: 16
    log:
        "logs/rule/merge_bams/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/merge_bams/{samples}_{build}_{method}.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        (echo "`date -R`: Merging BAM files..." &&
        samtools merge -f -@ {threads} {output.merged_bam} {input.bam1} {input.bam2} &&
        echo "`date -R`: Success! BAM files merged." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Sorting merged BAM..." &&
        samtools sort -@ {threads} -T {params.tmpdir} -o {output.merged_bam}.tmp {output.merged_bam} &&
        mv {output.merged_bam}.tmp {output.merged_bam} &&
        echo "`date -R`: Success! BAM sorted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Indexing merged BAM..." &&
        samtools index {output.merged_bam} {output.merged_bam_idx} &&
        echo "`date -R`: Success! BAM indexed." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """

rule bam2bed_se_filtered:
    input:
        rules.merge_bams.output.merged_bam,
    output:
        bed="results/{method}/{samples}/{samples}_{build}.bed",
        bam="results/{method}/{samples}/{samples}_{build}.bam",
        idx="results/{method}/{samples}/{samples}_{build}.bam.bai",
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
