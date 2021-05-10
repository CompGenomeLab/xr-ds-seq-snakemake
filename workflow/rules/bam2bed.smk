
rule bam2bed_se:
    input:
        "results/{samples}/{samples}_cutadapt_se_{build}.bam",
    output:
        "results/{samples}/{samples}_{build}_se.bed",
    params:
        q_trim=config["samtools_q_trim_se"], 
    log:
        "logs/{samples}/{samples}_{build}_bam2bed_se.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_bam2bed_se.benchmark.txt",
    conda:
        "../envs/bam2bed.yaml"
    shell:  
        """
        (echo "`date -R`: Processing bam file..." && 
        samtools view {params.q_trim} -b {input} |&
        bedtools bamtobed > {output} &&
        echo "`date -R`: Success! Bam file converted to bed format." || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        """

rule bam2bed_pe:
    input:
        "results/{samples}/{samples}_cutadapt_pe_{build}.bam",
    output:
        bed="results/{samples}/{samples}_{build}_pe.bed",
        bam=temp("results/{samples}/{samples}_{build}_sorted.bam"),
    params:
        q_trim=config["samtools_q_trim_pe"],
    log:
        "logs/{samples}/{samples}_{build}_bam2bed_pe.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_bam2bed_pe.benchmark.txt",
    conda:
        "../envs/bam2bed.yaml"
    shell:  
        """
        (echo "`date -R`: Sorting bam file..." &&
        samtools sort -n {input} > {output.bam} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

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
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """