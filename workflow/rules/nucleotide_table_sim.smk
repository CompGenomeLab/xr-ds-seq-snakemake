
rule nucleotide_table_ds_sim:
    input:
        "results/{samples}/{samples}_{build}_ds_sim.fa",
    output:
        dinuc="results/{samples}/{samples}_{build}_sim_dinucleotideTable.txt",
        nuc="results/{samples}/{samples}_{build}_sim_nucleotideTable.txt",
    log:
        "logs/{samples}/{samples}_{build}_nucleotide_table_ds_sim.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_nucleotide_table_ds_sim.benchmark.txt",
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

rule nucleotide_table_xr_sim:
    input:
        "results/{samples}/{samples}_{build}_xr_sim.fa",
    output:
        dinuc="results/{samples}/{samples}_{build}_sim_dinucleotideTable.txt",
        nuc="results/{samples}/{samples}_{build}_sim_nucleotideTable.txt",
        filt=temp("results/{samples}/{samples}_{build}_xr_sim_filt.fa"),
    log:
        "logs/{samples}/{samples}_{build}_nucleotide_table_xr_sim.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_nucleotide_table_xr_sim.benchmark.txt",
    shell:  
        """
        len="$(awk 'BEGIN{{RS=">";ORS=""}}{{print length($2)"\\n"}}' {input} | sort | uniq -c | sort -k1,1n | tail -1 | awk '{{print $2}}')"

        awk -v num="$len" 'BEGIN{{RS=">";ORS=""}}length($2)==num{{print ">"$0}}' {input} > {output.filt}

        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {output.filt} \
        -k 2 \
        -o {output.dinuc} 2>{log}

        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {output.filt} \
        -k 1 \
        -o {output.nuc} 2>>{log}
        """