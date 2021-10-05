
rule rename_raw_input:
    input:
        r1=lambda w: getPaired(w.samples, config["sample"], "forward", "resources/input/"),
        r2=lambda w: getPaired(w.samples, config["sample"], "reverse", "resources/input/"),
    output:
        r1="resources/input/{samples}_1.fastq.gz", 
        r2="resources/input/{samples}_2.fastq.gz",
    log:
        "logs/{samples}/{samples}_rename_raw_input.log",
    benchmark:
        "logs/{samples}/{samples}_rename_raw_input.benchmark.txt",
    shell:
        """
        mv {input.r1} {output.r1}

        mv {input.r2} {output.r2}
        """