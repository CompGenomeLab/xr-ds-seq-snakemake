
rule bam2bed_se:
    input:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_se.bam"
    output:
        "results/{dir}/{samples}{v}/{samples}_cutadapt.bed"
    params:
        q_trim=lambda w: getSampleInfo(w, config["samtools_q_trim"]) 
    log:
        "results/{dir}/{samples}{v}/log/bam2bed.log"
    benchmark:
        "results/{dir}/{samples}{v}/log/bam2bed.benchmark.txt",
    conda:
        "../envs/bam2bed.yaml"
    shell:  
        """
        (samtools view {params.q_trim} -b {input} |
        bedtools bamtobed > {output}) 2>{log}
        """


rule bam2bed_pe:
    input:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_pe.bam"
    output:
        "results/{dir}/{samples}{v}/{samples}_cutadapt.bed"
    params:
        q_trim=lambda w: getSampleInfo(w, config["samtools_q_trim"])
    log:
        "results/{dir}/{samples}{v}/log/bam2bed.log"
    benchmark:
        "results/{dir}/{samples}{v}/log/bam2bed.benchmark.txt",
    conda:
        "../envs/bam2bed.yaml"
    shell:  
        """
        samtools sort -n {input} |
        samtools view {params.q_trim} - |
        bedtools bamtobed -bedpe -mate1 |
        awk '{{\
            if ($9=="+")\
                print $1"\t"$2"\t"$6"\t"$7"\t"$8"\t"$9;\
            else if ($9=="-")\
                print $1"\t"$5"\t"$3"\t"$7"\t"$8"\t"$9;\
            }}' > {output}
        """

