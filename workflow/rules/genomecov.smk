
rule genomecov_ds:
    input:
        plus="results/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        minus="results/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed",
        ref_genome="resources/ref_genomes/{build}/genome_{build}.fa.fai",
    output:
        plus=temp("results/{samples}/{samples}_{build}_DS_sorted_plus.bdg"),
        minus=temp("results/{samples}/{samples}_{build}_DS_sorted_minus.bdg"),
    params:
        read=lambda w, input: mappedReads(input),
    log:
        "logs/{samples}/{samples}_{build}_genomecov_ds.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_genomecov_ds.benchmark.txt",
    conda:
        "../envs/genomecov.yaml"
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

rule genomecov_xr:
    input:
        plus=temp("results/{samples}/{samples}_{build}_sorted_plus.bed"),
        minus=temp("results/{samples}/{samples}_{build}_sorted_minus.bed"),
        ref_genome="resources/ref_genomes/{build}/genome_{build}.fa.fai",
    output:
        plus="results/{samples}/{samples}_{build}_XR_sorted_plus.bdg",
        minus="results/{samples}/{samples}_{build}_XR_sorted_minus.bdg",
    params:
        read=lambda w, input: mappedReads(input),
    log:
        "logs/{samples}/{samples}_{build}_genomecov_xr.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_genomecov_xr.benchmark.txt",
    conda:
        "../envs/genomecov.yaml"
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