
rule genomecov_ds:
    input:
        plus="results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        minus="results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed",
        ref_genome="resources/ref_genomes/{build}/genome_{build}.fa.fai",
    output:
        plus=temp("results/{method}/{samples}/{samples}_{build}_DS_sorted_plus.bdg"),
        minus=temp("results/{method}/{samples}/{samples}_{build}_DS_sorted_minus.bdg"),
    params:
        read=lambda w, input: mappedReads(input[0], input[1]),
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_genomecov_ds.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_genomecov_ds.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Calculating genome coverage of {input.plus}..." &&
        bedtools genomecov \
        -i {input.plus} \
        -g {input.ref_genome} \
        -bg \
        -scale $(echo {params.read} | awk '{{print 1000000/$1}}') \
        > {output.plus} &&
        echo "`date -R`: Success! Genome coverage is calculated." || 
        {{ echo "`date -R`: Process failed..."; rm {output.plus}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Calculating genome coverage of {input.minus}..." &&
        bedtools genomecov \
        -i {input.minus} \
        -g {input.ref_genome} \
        -bg \
        -scale $(echo {params.read} | awk '{{print 1000000/$1}}') \
        > {output.minus} &&
        echo "`date -R`: Success! Genome coverage is calculated." || 
        {{ echo "`date -R`: Process failed..."; rm {output.minus}; exit 1; }}  ) >> {log} 2>&1
        """

rule genomecov_xr:
    input:
        plus="results/{method}/{samples}/{samples}_{build}_sorted_plus.bed",
        minus="results/{method}/{samples}/{samples}_{build}_sorted_minus.bed",
        ref_genome="resources/ref_genomes/{build}/genome_{build}.fa.fai",
    output:
        plus=temp("results/{method}/{samples}/{samples}_{build}_XR_sorted_plus.bdg"),
        minus=temp("results/{method}/{samples}/{samples}_{build}_XR_sorted_minus.bdg"),
    params:
        read=lambda w, input: mappedReads(input[0], input[1]),
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_genomecov_xr.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_genomecov_xr.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Calculating genome coverage of {input.plus}..." &&
        bedtools genomecov \
        -i {input.plus} \
        -g {input.ref_genome} \
        -bg \
        -scale $(echo {params.read} | awk '{{print 1000000/$1}}') \
        > {output.plus} &&
        echo "`date -R`: Success! Genome coverage is calculated." || 
        echo "`date -R`: Process failed...") > {log} 2>&1


        (echo "`date -R`: Calculating genome coverage of {input.minus}..." &&
        bedtools genomecov \
        -i {input.minus} \
        -g {input.ref_genome} \
        -bg \
        -scale $(echo {params.read} | awk '{{print 1000000/$1}}') \
        > {output.minus} &&
        echo "`date -R`: Success! Genome coverage is calculated." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """