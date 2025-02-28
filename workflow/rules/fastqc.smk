
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