
rule bedGraphToBigWig:
    input:
        bdg="results/{method}/{samples}/{samples}_{build}_{method}_sorted_{strand}.bdg",
        index="resources/ref_genomes/{build}/genome_{build}.fa.fai",
    output:
        sort=temp("results/{method}/{samples}/{samples}_{build}_{method}_resorted_{strand}.bdg"),
        bw=report("results/{method}/{samples}/{samples}_{build}_{method}_sorted_{strand}.bw", 
                category="BigWig"),
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_bedGraphToBigWig_{strand}.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_bedGraphToBigWig_{strand}.benchmark.txt",
    conda:
        "../envs/bedGraphToBigWig.yaml"
    shell:  
        """
        (echo "`date -R`: Converting bedGraph to bigWig..." &&
        LC_COLLATE=C sort -k1,1 -k2,2n {input.bdg} > {output.sort} &&
        bedGraphToBigWig {output.sort} {input.index} {output.bw} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
        """