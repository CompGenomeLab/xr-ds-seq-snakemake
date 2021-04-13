
rule nucleotide_table_ds:
    input:
        "results/{samples}/{samples}_{build}_sorted_10.fa",
    output:
        dinuc=temp("results/{samples}/{samples}_{build}_sorted_dinucleotideTable.txt"),
        nuc=temp("results/{samples}/{samples}_{build}_sorted_nucleotideTable.txt"),
    log:
        "logs/{samples}/{samples}_{build}_nucleotide_table_ds.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_nucleotide_table_ds.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Calculating dinucleotide abundance table..." &&
        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {input} \
        -k 2 \
        -o {output.dinuc} &&
        echo "`date -R`: Success! Dinucleotide abundance table is calculated." ||
        echo "`date -R`: Dinucleotide abundace table cannot be calculated...") \
        > {log} 2>&1

        (echo "`date -R`: Calculating nucleotide abundance table..." &&
        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {input} \
        -k 1 \
        -o {output.nuc}  &&
        echo "`date -R`: Success! Nucleotide abundance table is calculated." ||
        echo "`date -R`: Nucleotide abundace table cannot be calculated...") \
        >> {log} 2>&1
        """

rule nucleotide_table_xr:
    input:
        "results/{samples}/{samples}_{build}_lengthMode.fa",
    output:
        dinuc=temp("results/{samples}/{samples}_{build}_sorted_dinucleotideTable.txt"),
        nuc=temp("results/{samples}/{samples}_{build}_sorted_nucleotideTable.txt"),
    log:
        "logs/{samples}/{samples}_{build}_nucleotide_table_xr.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_nucleotide_table_xr.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Calculating dinucleotide abundance table..." &&
        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {input} \
        -k 2 \
        -o {output.dinuc}  &&
        echo "`date -R`: Success! Dinucleotide abundance table is calculated." ||
        echo "`date -R`: Dinucleotide abundace table cannot be calculated...") \
        > {log} 2>&1

        (echo "`date -R`: Calculating nucleotide abundance table..." &&
        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {input} \
        -k 1 \
        -o {output.nuc}  &&
        echo "`date -R`: Success! Nucleotide abundance table is calculated." ||
        echo "`date -R`: Nucleotide abundace table cannot be calculated...") \
        >> {log} 2>&1
        """