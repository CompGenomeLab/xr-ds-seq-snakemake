
rule genomecov_ds:
    input:
        bed="results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_{strand}.bed",
        ref_genome=rules.genome_indexing.output,
    output:
        temp("results/{method}/{samples}/{samples}_{build}_DS_sorted_{strand}.bdg"),
    params:
        read=lambda w, input: mappedReads(input[0]),
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_genomecov_ds_{strand}.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_genomecov_ds_{strand}.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Calculating genome coverage of {input.bed}..." &&
        bedtools genomecov \
        -i {input.bed} \
        -g {input.ref_genome} \
        -bg \
        -scale $(echo {params.read} | awk '{{print 1000000/$1}}') \
        > {output} &&
        echo "`date -R`: Success! Genome coverage is calculated." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
        """

rule genomecov_xr:
    input:
        bed="results/{method}/{samples}/{samples}_{build}_sorted_{strand}.bed",
        ref_genome=rules.genome_indexing.output,
    output:
        bed=temp("results/{method}/{samples}/{samples}_{build}_XR_sorted_{strand}.bdg"),
    params:
        read=lambda w, input: mappedReads(input[0]),
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_genomecov_xr_{strand}.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_genomecov_xr_{strand}.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Calculating genome coverage of {input.bed}..." &&
        bedtools genomecov \
        -i {input.bed} \
        -g {input.ref_genome} \
        -bg \
        -scale $(echo {params.read} | awk '{{print 1000000/$1}}') \
        > {output.bed} &&
        echo "`date -R`: Success! Genome coverage is calculated." || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        """