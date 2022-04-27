
rule reposition:
    input:
        plus=rules.sep_strands.output.plus,
        minus=rules.sep_strands.output.minus,
        index=rules.genome_indexing.output,
    output:
        plus=temp("results/{method}/{samples}/{samples}_{build}_sorted_plus_10.bed"),
        minus=temp("results/{method}/{samples}/{samples}_{build}_sorted_minus_10.bed"), 
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_reposition.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_reposition.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Centering {input.plus} at damage site..." &&
        bedtools flank \
        -i {input.plus} \
        -g {input.index} \
        -l 6 \
        -r 0 |& 
        bedtools slop \
        -g {input.index} \
        -l 0 \
        -r 4 |& 
        awk '{{ if ($3-$2 == 10) {{ print }} }}' \
        > {output.plus} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.plus}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Centering {input.minus} at damage site..." &&
        bedtools flank \
        -i {input.minus} \
        -g {input.index} \
        -l 0 \
        -r 6 |& 
        bedtools slop \
        -g {input.index} \
        -l 4 \
        -r 0 |& 
        awk '{{ if ($3-$2 == 10) {{ print }} }}' \
        > {output.minus} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.minus}; exit 1; }}  ) >> {log} 2>&1
        """