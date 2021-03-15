
rule dinucleotide_table_ds:
    input:
        "results/{samples}/{samples}_{build}_sorted_10.fa",
    output:
        "results/{samples}/{samples}_{build}_sorted_dinucleotideTable.txt",
    params:
        k="2",
    log:
        "logs/{samples}/{samples}_{build}_dinucleotide_table_ds.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_dinucleotide_table_ds.benchmark.txt",
    shell:  
        """
        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {input} \
        -k {params.k} \
        --percentage \
        -o {output} 2>{log}
        """

rule dinucleotide_table_xr:
    input:
        "results/{samples}/{samples}_{build}_lengthMode.fa",
    output:
        "results/{samples}/{samples}_{build}_sorted_dinucleotideTable.txt",
    params:
        k="2",
    log:
        "logs/{samples}/{samples}_{build}_dinucleotide_table_xr.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_dinucleotide_table_xr.benchmark.txt",
    shell:  
        """
        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {input} \
        -k {params.k} \
        -o {output} 2>{log}
        """

