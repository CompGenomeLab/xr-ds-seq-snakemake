
rule reposition:
    input:
        plus="results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_plus.bed",
        minus="results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_minus.bed",
    output:
        plus="results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_plus_10.bed",
        minus="results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_minus_10.bed", 
    params:
        index=lambda w: getGenome(w, "index"),    
    log:
        "results/{dir}/{samples}{v}/log/reposition.log",
    benchmark:
        "results/{dir}/{samples}{v}/log/reposition.benchmark.txt",
    conda:
        "../envs/reposition.yaml"
    shell:  
        """
        bedtools flank -i {input.plus} -g {params.index} -l 6 -r 0 | 
        bedtools slop -g {params.index} -l 0 -r 4 | 
        awk '{{ if ($3-$2 == 10) {{ print }} }}' > {output.plus}

        bedtools flank -i {input.minus} -g {params.index} -l 0 -r 6 | 
        bedtools slop -g {params.index} -l 4 -r 0 | 
        awk '{{ if ($3-$2 == 10) {{ print }} }}' > {output.minus}
        """

