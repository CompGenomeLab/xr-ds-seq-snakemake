
rule dinucleotide_table_ds:
    input:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_10.fa",
    output:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_dinucleotideTable.txt",
    params:
        k="2",
        l="",
        percentage="",
    log:
        "results/{dir}/{samples}{v}/log/dinucleotide_table_ds.log",
    benchmark:
        "results/{dir}/{samples}{v}/log/dinucleotide_table_ds.benchmark.txt",
    script:  
        "../scripts/fa2kmerAbundanceTable.py" 

rule dinucleotide_table_xr:
    input:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_lengthMode.fa",
    output:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_dinucleotideTable.txt",
    params:
        k="2",
        l="",
        percentage="",
    log:
        "results/{dir}/{samples}{v}/log/dinucleotide_table_xr.log",
    benchmark:
        "results/{dir}/{samples}{v}/log/dinucleotide_table_xr.benchmark.txt",
    script:  
        "../scripts/fa2kmerAbundanceTable.py" 

