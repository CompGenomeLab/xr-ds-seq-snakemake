
rule zip_samples:
    input:
        "resources/samples/{samples}",
    output:
        "resources/samples/{samples}.gz",
    log:
        "logs/{samples}/{samples}_zip_samples.log",
    benchmark:
        "logs/{samples}/{samples}_zip_samples.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Zipping raw data..." && 
        gzip {input} &&
        echo "`date -R`: Success! Zipping is done." ||
        echo "`date -R`: Zipping failed...") > {log} 2>&1
        """