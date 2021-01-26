
rule genomecov_ds:
    input:
        plus="results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_plus_dipyrimidines.bed",
        minus="results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_minus_dipyrimidines.bed",
    output:
        plus="results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_plus_dipyrimidines.bdg",
        minus="results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_minus_dipyrimidines.bdg",
    params:
        read=lambda w, input: mappedReads(input),
        ref_genome=lambda w: getGenome(w, "index"),
    log:
        "results/{dir}/{samples}{v}/log/genomecov_ds.log",
    benchmark:
        "results/{dir}/{samples}{v}/log/genomecov_ds.benchmark.txt",
    conda:
        "../envs/genomecov.yaml"
    shell:  
        """
        bedtools genomecov -i {input.plus} -g {params.ref_genome} -bg -scale \
        $(echo {params.read} | awk '{{print 1000000/$1}}') > {output.plus}

        bedtools genomecov -i {input.minus} -g {params.ref_genome} -bg -scale \
        $(echo {params.read} | awk '{{print 1000000/$1}}') > {output.minus}
        """

rule genomecov_xr:
    input:
        plus="results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_plus.bed",
        minus="results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_minus.bed",
    output:
        plus="results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_plus.bdg",
        minus="results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_minus.bdg",
    params:
        read=lambda w, input: mappedReads(input),
        ref_genome=lambda w: getGenome(w, "index"),
    log:
        "results/{dir}/{samples}{v}/log/genomecov_xr.log",
    benchmark:
        "results/{dir}/{samples}{v}/log/genomecov_xr.benchmark.txt",
    conda:
        "../envs/genomecov.yaml"
    shell:  
        """
        bedtools genomecov -i {input.plus} -g {params.ref_genome} -bg -scale \
        $(echo {params.read} | awk '{{print 1000000/$1}}') > {output.plus}

        bedtools genomecov -i {input.minus} -g {params.ref_genome} -bg -scale \
        $(echo {params.read} | awk '{{print 1000000/$1}}') > {output.minus}
        """