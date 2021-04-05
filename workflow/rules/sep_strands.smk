
rule sep_strands:
    input:
        "results/{samples}/{samples}_{build}_sorted_chr.bed",
    output:
        plus="results/{samples}/{samples}_{build}_sorted_plus.bed",
        minus="results/{samples}/{samples}_{build}_sorted_minus.bed",
    log:
        "logs/{samples}/{samples}_{build}_sep_strands.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_sep_strands.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Separating plus stranded reads..." &&
        awk '{{if($6=="+"){{print}}}}' {input} > {output.plus} &&
        echo "`date -R`: Success! Reads are separated." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Separating minus stranded reads..." &&
        awk '{{if($6=="-"){{print}}}}' {input} > {output.minus} &&
        echo "`date -R`: Success! Reads are separated." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """