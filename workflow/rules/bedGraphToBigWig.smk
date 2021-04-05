
rule bedGraphToBigWig_ds:
    input:
        bdg="results/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_{strand}.bdg",
        index="resources/ref_genomes/{build}/genome_{build}.fa.fai",
    output:
        "results/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_{strand}.bw",
    log:
        "logs/{samples}/{samples}_{build}_bedGraphToBigWig_{strand}.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_bedGraphToBigWig_{strand}.benchmark.txt",
    conda:
        "../envs/bedGraphToBigWig.yaml"
    shell:  
        """
        (echo "`date -R`: Converting bedGraph to bigWig..." &&
        bedGraphToBigWig {input.bdg} {input.index} {output} &&
        echo "`date -R`: Success! Conversion is done." || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        """

rule bedGraphToBigWig_xr:
    input:
        bdg="results/{samples}/{samples}_{build}_sorted_xr_{strand}.bdg",
        index="resources/ref_genomes/{build}/genome_{build}.fa.fai", 
    output:
        "results/{samples}/{samples}_{build}_sorted_xr_{strand}.bw",
    log:
        "logs/{samples}/{samples}_{build}_bedGraphToBigWig_{strand}.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_bedGraphToBigWig_{strand}.benchmark.txt",
    conda:
        "../envs/bedGraphToBigWig.yaml"
    shell:  
        """
        (echo "`date -R`: Converting bedGraph to bigWig..." &&
        bedGraphToBigWig {input.bdg} {input.index} {output} &&
        echo "`date -R`: Success! Conversion is done." || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        """