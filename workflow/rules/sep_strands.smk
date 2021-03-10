
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
        awk '{{if($6=="+"){{print}}}}' {input} > {output.plus}
        awk '{{if($6=="-"){{print}}}}' {input} > {output.minus}
        """
