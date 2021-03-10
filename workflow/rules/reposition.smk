
rule reposition:
    input:
        plus="results/{samples}/{samples}_{build}_sorted_plus.bed",
        minus="results/{samples}/{samples}_{build}_sorted_minus.bed",
        index="resources/ref_genomes/{build}/genome_{build}.fa.fai",
    output:
        plus="results/{samples}/{samples}_{build}_sorted_plus_10.bed",
        minus="results/{samples}/{samples}_{build}_sorted_minus_10.bed", 
    log:
        "logs/{samples}/{samples}_{build}_reposition.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_reposition.benchmark.txt",
    conda:
        "../envs/reposition.yaml"
    shell:  
        """
        bedtools flank -i {input.plus} -g {input.index} -l 6 -r 0 | 
        bedtools slop -g {input.index} -l 0 -r 4 | 
        awk '{{ if ($3-$2 == 10) {{ print }} }}' > {output.plus}

        bedtools flank -i {input.minus} -g {input.index} -l 0 -r 6 | 
        bedtools slop -g {input.index} -l 4 -r 0 | 
        awk '{{ if ($3-$2 == 10) {{ print }} }}' > {output.minus}
        """

