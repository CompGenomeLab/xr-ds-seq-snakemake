rule rename_raw:
    input:
        r1="resources/samples/{samples}_R1.fastq.gz", 
        r2="resources/samples/{samples}_R2.fastq.gz",
    output:
        r1=temp("resources/samples/{samples}_1.fastq.gz"), 
        r2=temp("resources/samples/{samples}_2.fastq.gz"),
    log:
        "logs/{samples}/{samples}_rename_raw.log",
    benchmark:
        "logs/{samples}/{samples}_rename_raw.benchmark.txt",
    shell:
        """
        cp {input.r1} {output.r1}

        cp {input.r2} {output.r2}
        """
    
