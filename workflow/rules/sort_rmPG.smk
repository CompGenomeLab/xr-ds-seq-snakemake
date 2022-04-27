rule sort_rmPG_se:
    input:
        rules.bowtie2_se.output.bam,
    output:
        header=temp("results/{method}/{samples}/{samples}_cutadapt_se_{build}_header.txt"),
        sort=temp("results/{method}/{samples}/{samples}_cutadapt_se_{build}_samSorted.bam"),
    log:
        "logs/rule/analysis/{samples}/{samples}_{method}_{build}_se_sort_rmPG.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{method}_{build}_se_sort_rmPG.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        (echo "`date -R`: Parsing the headers..." &&
        samtools view -H {input} | grep -v "^@PG" > {output.header} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.header}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Parsing the headers..." &&
        samtools reheader {output.header} {input} | samtools sort -o {output.sort} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.sort}; exit 1; }}  ) >> {log} 2>&1
        """

rule sort_rmPG_pe:
    input:
        rules.bowtie2_pe.output.bam,
    output:
        header=temp("results/{method}/{samples}/{samples}_cutadapt_pe_{build}_header.txt"),
        sort=temp("results/{method}/{samples}/{samples}_cutadapt_pe_{build}_samSorted.bam"),
    log:
        "logs/rule/analysis/{samples}/{samples}_{method}_{build}_pe_sort_rmPG.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{method}_{build}_pe_sort_rmPG.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        (echo "`date -R`: Parsing the headers excluding the @PG headers..." &&
        samtools view -H {input} | grep -v "^@PG" > {output.header} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.header}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Reheader and sort..." &&
        samtools reheader {output.header} {input} | samtools sort -o {output.sort} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.sort}; exit 1; }}  ) >> {log} 2>&1
        """