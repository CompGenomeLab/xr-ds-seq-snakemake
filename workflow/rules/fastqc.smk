
rule fastqc_se:
    input:
        "resources/samples/{samples}.fastq.gz",
    output:
        html=report("results/processed_files/{samples}_fastqc.html", category="QC"),
        zip="results/processed_files/{samples}_fastqc.zip",
    wildcard_constraints:
        samples='|'.join([s for s in config["meta"].keys()]),
    params: 
        extra="",
        tmpdir="results/processed_files/{samples}",
    log:
        "logs/rule/fastqc_se/{samples}.log",
    benchmark:
        "logs/rule/fastqc_se/{samples}.benchmark.txt",
    threads: 16
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        (echo "`date -R`: FastQC..." &&
        mkdir -p {params.tmpdir} &&
        fastqc \
        {params.extra} \
        -t {threads} \
        --outdir {params.tmpdir} \
        {input} &&
        mv {params.tmpdir}/{wildcards.samples}_fastqc.html {output.html} && 
        mv {params.tmpdir}/{wildcards.samples}_fastqc.zip {output.zip} &&
        rmdir {params.tmpdir} && 
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1 
        """

rule fastqc_se_after_cutadapt:
    input:
        "results/XR/{samples}/{samples}_cutadapt.fastq.gz",
    output:
        html=report("results/processed_files/{samples}_cutadapt_fastqc.html", category="QC"),
        zip="results/processed_files/{samples}_cutadapt_fastqc.zip",
    wildcard_constraints:
        samples='|'.join([s for s in config["meta"].keys()]),
    params: 
        extra="",
        tmpdir="results/processed_files",
    log:
        "logs/rule/fastqc_se_after_cutadapt/{samples}_cutadapt.log",
    benchmark:
        "logs/rule/fastqc_se_after_cutadapt/{samples}_cutadapt.benchmark.txt",
    threads: 16
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        (echo "`date -R`: FastQC..." &&
        fastqc \
        {params.extra} \
        -t {threads} \
        --outdir {params.tmpdir} \
        {input} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1 
        """

rule fastqc_pe:
    input:
        "resources/samples/{samples}_{ext}.fastq.gz", 
    output:
        html=report("results/processed_files/{samples}_{ext}_fastqc.html", category="QC"), 
        zip="results/processed_files/{samples}_{ext}_fastqc.zip", 
    wildcard_constraints:
        samples='|'.join([s for s in config["meta"].keys()]),
        ext='R1|R2|1|2',
    params: 
        extra="",
        tmpdir="results/processed_files/{samples}_{ext}",
        name="{samples}_{ext}",
    log:
        "logs/rule/fastqc_pe/{samples}_{ext}.log",
    benchmark:
        "logs/rule/fastqc_pe/{samples}_{ext}.benchmark.txt",
    threads: 16
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        (echo "`date -R`: FastQC..." &&
        mkdir -p {params.tmpdir} &&
        fastqc \
        {params.extra} \
        -t {threads} \
        --outdir {params.tmpdir} \
        {input} &&
        mv {params.tmpdir}/{params.name}_fastqc.html {output.html} && 
        mv {params.tmpdir}/{params.name}_fastqc.zip {output.zip} &&
        rmdir {params.tmpdir} && 
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1 
        """

rule fastqc_unmapped:
    input:
        "results/{method}/{samples}/{samples}_cutadapt_se_{build}_unmapped.fastq"
    output:
        html=report("results/{method}/{samples}/{samples}_cutadapt_se_{build}_unmapped_fastqc.html", category="QC"),
        zip="results/{method}/{samples}/{samples}_cutadapt_se_{build}_unmapped_fastqc.zip"
    log:
        "logs/rule/fastqc_unmapped/{samples}_{build}_{method}.log"
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        (echo "`date -R`: Running FastQC on unmapped reads..." &&
        fastqc {input} -o $(dirname {output.html}) -t {threads} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """

rule fastqc_unmapped_filtered:
    input:
        "results/{method}/{samples}/{samples}_cutadapt_se_{build}_unmapped_filtered_overrep.fastq"
    output:
        html=report("results/{method}/{samples}/{samples}_cutadapt_se_{build}_unmapped_filtered_overrep_fastqc.html", category="QC"),
        zip="results/{method}/{samples}/{samples}_cutadapt_se_{build}_unmapped_filtered_overrep_fastqc.zip"
    log:
        "logs/rule/fastqc_unmapped_filtered/{samples}_{build}_{method}.log"
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        (echo "`date -R`: Running FastQC on unmapped reads..." &&
        fastqc {input} -o $(dirname {output.html}) -t {threads} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """