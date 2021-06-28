rule rename_raw:
    input:
        r1=lambda w: getPaired(w.samples, config["sample"], "forward"),
        r2=lambda w: getPaired(w.samples, config["sample"], "reverse"),
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
    
