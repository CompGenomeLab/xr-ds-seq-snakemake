rule rename_raw:
    input:
        r1=lambda w: getPaired(w.samples, config["sample"], "forward", "resources/samples/"),
        r2=lambda w: getPaired(w.samples, config["sample"], "reverse", "resources/samples/"),
    output:
        r1="resources/samples/{samples}_1.fastq.gz", 
        r2="resources/samples/{samples}_2.fastq.gz",
    log:
        "logs/rule/analysis/{samples}/{samples}_rename_raw.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_rename_raw.benchmark.txt",
    shell:
        """
        mv {input.r1} {output.r1}

        mv {input.r2} {output.r2}
        """