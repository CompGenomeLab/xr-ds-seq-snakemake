
rule nucleotide_table:
    input:
        lambda w: input4nucTable(w, config["meta"]),
    output:
        dinuc=temp("results/{method}/{samples}/{samples}_{build}_sorted_dinucleotideTable.txt"),
        nuc=temp("results/{method}/{samples}/{samples}_{build}_sorted_nucleotideTable.txt"),
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_nucleotide_table.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_nucleotide_table.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Calculating dinucleotide abundance table..." &&
        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {input} \
        -k 2 \
        -o {output.dinuc} &&
        echo "`date -R`: Success! Dinucleotide abundance table is calculated." ||
        {{ echo "`date -R`: Dinucleotide abundace table cannot be calculated..."; rm {output.dinuc}; exit 1; }}  ) \
        > {log} 2>&1

        (echo "`date -R`: Calculating nucleotide abundance table..." &&
        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {input} \
        -k 1 \
        -o {output.nuc}  &&
        echo "`date -R`: Success! Nucleotide abundance table is calculated." ||
        {{ echo "`date -R`: Nucleotide abundace table cannot be calculated..."; rm {output.nuc}; exit 1; }}  ) \
        >> {log} 2>&1
        """