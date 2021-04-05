
rule nucleotide_table_ds:
    input:
        "results/{samples}/{samples}_{build}_sorted_10.fa",
    output:
        dinuc="results/{samples}/{samples}_{build}_sorted_dinucleotideTable.txt",
        nuc="results/{samples}/{samples}_{build}_sorted_nucleotideTable.txt",
    log:
        "logs/{samples}/{samples}_{build}_nucleotide_table_ds.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_nucleotide_table_ds.benchmark.txt",
    shell:  
        """
        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {input} \
        -k 2 \
        -o {output.dinuc} 2>{log}

        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {input} \
        -k 1 \
        -o {output.nuc} 2>>{log}
        """

rule nucleotide_table_xr:
    input:
        "results/{samples}/{samples}_{build}_lengthMode.fa",
    output:
        dinuc="results/{samples}/{samples}_{build}_sorted_dinucleotideTable.txt",
        nuc="results/{samples}/{samples}_{build}_sorted_nucleotideTable.txt",
    log:
        "logs/{samples}/{samples}_{build}_nucleotide_table_xr.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_nucleotide_table_xr.benchmark.txt",
    shell:  
        """
        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {input} \
        -k 2 \
        -o {output.dinuc} 2>{log}

        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {input} \
        -k 1 \
        -o {output.nuc} 2>>{log}
        """