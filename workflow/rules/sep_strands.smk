
rule sep_strands:
    input:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_chr.bed",
    output:
        plus="results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_plus.bed",
        minus="results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_minus.bed",
    log:
        "results/{dir}/{samples}{v}/log/sep_strands.log",
    benchmark:
        "results/{dir}/{samples}{v}/log/sep_strands.benchmark.txt",
    shell:  
        """
        awk '{{if($6=="+"){{print}}}}' {input} > {output.plus}
        awk '{{if($6=="-"){{print}}}}' {input} > {output.minus}
        """
