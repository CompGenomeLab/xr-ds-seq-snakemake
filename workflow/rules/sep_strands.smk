
rule sep_strands:
    input:
       rules.sort_filter.output,
    output:
        plus="results/{method}/{samples}/{samples}_{build}_sorted_plus.bed",
        minus="results/{method}/{samples}/{samples}_{build}_sorted_minus.bed",
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_sep_strands.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_sep_strands.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Separating plus stranded reads..." &&
        awk '{{if($6=="+"){{print}}}}' {input} > {output.plus} &&
        echo "`date -R`: Success! Reads are separated." || 
        {{ echo "`date -R`: Process failed..."; rm {output.plus}; exit 1; }}  )  > {log} 2>&1

        (echo "`date -R`: Separating minus stranded reads..." &&
        awk '{{if($6=="-"){{print}}}}' {input} > {output.minus} &&
        echo "`date -R`: Success! Reads are separated." || 
        {{ echo "`date -R`: Process failed..."; rm {output.minus}; exit 1; }}  )  >> {log} 2>&1
        """