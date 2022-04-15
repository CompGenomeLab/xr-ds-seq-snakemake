
rule bed2fasta_ds:
    input:
        plus="results/{method}/{samples}/{samples}_{build}_sorted_plus_10.bed",
        minus="results/{method}/{samples}/{samples}_{build}_sorted_minus_10.bed",
        genome="resources/ref_genomes/{build}/genome_{build}.fa",
    output:
        plus=temp("results/{method}/{samples}/{samples}_{build}_sorted_plus_10.fa"),
        minus=temp("results/{method}/{samples}/{samples}_{build}_sorted_minus_10.fa"),
        comb=temp("results/{method}/{samples}/{samples}_{build}_sorted_10.fa"),
        bed=temp("results/{method}/{samples}/{samples}_{build}_sorted_10.bed"),       
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_bed2fasta_ds.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_bed2fasta_ds.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        (echo "`date -R`: Combine files..." &&
        cat {input.plus} {input.minus} > {output.bed} &&
        echo "`date -R`: Success! Files are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        
        (echo "`date -R`: Converting {output.bed} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {output.bed} \
        -fo {output.comb} \
        -s &&
        echo "`date -R`: Success! {output.bed} is converted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Converting {input.plus} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.plus} \
        -fo {output.plus} \
        -s &&
        echo "`date -R`: Success! {input.plus} is converted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Converting {input.minus} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.minus} \
        -fo {output.minus} \
        -s &&
        echo "`date -R`: Success! {input.minus} is converted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """

rule bed2fasta_xr:
    input:
        bed="results/{method}/{samples}/{samples}_{build}_lengthMode.bed",
        genome="resources/ref_genomes/{build}/genome_{build}.fa",
    output:
        temp("results/{method}/{samples}/{samples}_{build}_lengthMode.fa"),
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_bed2fasta_xr.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_bed2fasta_xr.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        (echo "`date -R`: Converting {input.bed} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.bed} \
        -fo {output} \
        -s &&
        echo "`date -R`: Success! {input.bed} is converted." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
        """

rule bed2fasta_input:
    input:
        bed=lambda w: input4inpFasta(w, config["meta"], config["sample"]),
        genome="resources/ref_genomes/{build}/genome_{build}.fa",
    output:
        "results/input/{samples}/{samples}_{build}.fasta",
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_bed2fasta_input.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_bed2fasta_input.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        (echo "`date -R`: Converting {input.bed} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.bed} \
        -fo {output} \
        -s &&
        echo "`date -R`: Success! {input.bed} is converted." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
        """