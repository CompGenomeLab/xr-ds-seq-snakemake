
rule bed2fasta_ds:
    input:
        plus="results/{samples}/{samples}_{build}_sorted_plus_10.bed",
        minus="results/{samples}/{samples}_{build}_sorted_minus_10.bed",
        genome="resources/ref_genomes/{build}/genome_{build}.fa",
    output:
        plus="results/{samples}/{samples}_{build}_sorted_plus_10.fa",
        minus="results/{samples}/{samples}_{build}_sorted_minus_10.fa",
        comb="results/{samples}/{samples}_{build}_sorted_10.fa",
        bed=temp("results/{samples}/{samples}_{build}_sorted_10.bed"),       
    log:
        "logs/{samples}/{samples}_{build}_bed2fasta_ds.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_bed2fasta_ds.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        (echo "`date -R`: Combine files..." &&
        cat {input.plus} {input.minus} > {output.bed} &&
        echo "`date -R`: Success! Files are combined." || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        
        (echo "`date -R`: Converting {output.bed} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {output.bed} \
        -fo {output.comb} \
        -s &&
        echo "`date -R`: Success! {output.bed} is converted." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Converting {input.plus} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.plus} \
        -fo {output.plus} \
        -s &&
        echo "`date -R`: Success! {input.plus} is converted." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Converting {input.minus} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.minus} \
        -fo {output.minus} \
        -s &&
        echo "`date -R`: Success! {input.minus} is converted." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule bed2fasta_xr:
    input:
        bed="results/{samples}/{samples}_{build}_lengthMode.bed",
        genome="resources/ref_genomes/{build}/genome_{build}.fa",
    output:
        "results/{samples}/{samples}_{build}_lengthMode.fa",
    log:
        "logs/{samples}/{samples}_{build}_bed2fasta_xr.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_bed2fasta_xr.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        (echo "`date -R`: Converting {input.bed} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.bed} \
        -fo {output} \
        -s &&
        echo "`date -R`: Success! {input.bed} is converted." || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        """