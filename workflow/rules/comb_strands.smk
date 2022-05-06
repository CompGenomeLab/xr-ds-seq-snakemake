
rule comb_strands:
    input:
        plus="results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        minus="results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed",
    output:
        "results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines.bed",
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_comb_strands.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_comb_strands.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Combine files..." &&
        cat {input.plus} {input.minus} > {output} &&
        echo "`date -R`: Success! Files are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """