
rule nucleotide_table_ds_sim:
    input:
        "results/{samples}/{samples}_{build}_ds_sim.fa",
    output:
        dinuc=temp("results/{samples}/{samples}_{build}_DS_sim_dinucleotideTable.txt"),
        nuc=temp("results/{samples}/{samples}_{build}_DS_sim_nucleotideTable.txt"),
    log:
        "logs/{samples}/{samples}_{build}_nucleotide_table_ds_sim.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_nucleotide_table_ds_sim.benchmark.txt",
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
        -o {output.nuc} &&
        echo "`date -R`: Success! Nucleotide abundance table is calculated." ||
        echo "`date -R`: Nucleotide abundace table cannot be calculated...") \
        >> {log} 2>&1
        """

rule nucleotide_table_xr_sim:
    input:
        "results/{samples}/{samples}_{build}_xr_sim.fa",
    output:
        dinuc=temp("results/{samples}/{samples}_{build}_XR_sim_dinucleotideTable.txt"),
        nuc=temp("results/{samples}/{samples}_{build}_XR_sim_nucleotideTable.txt"),
        filt=temp("results/{samples}/{samples}_{build}_XR_sim_filt.fa"),
    log:
        "logs/{samples}/{samples}_{build}_nucleotide_table_xr_sim.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_nucleotide_table_xr_sim.benchmark.txt",
    shell:  
        """
        len="$(awk 'BEGIN{{RS=">";ORS=""}}{{print length($2)"\\n"}}' {input} | sort | uniq -c | sort -k1,1n | tail -1 | awk '{{print $2}}')"

        (echo "`date -R`: Filtering by most occurred read length..." &&
        awk -v num="$len" \
        'BEGIN{{RS=">";ORS=""}}length($2)==num{{print ">"$0}}' \
        {input} > {output.filt} &&
        echo "`date -R`: Success!" ||
        echo "`date -R`: Filtering is failed...") \
        > {log} 2>&1
        
        (echo "`date -R`: Calculating dinucleotide abundance table..." &&
        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {output.filt} \
        -k 2 \
        -o {output.dinuc} &&
        echo "`date -R`: Success! Dinucleotide abundance table is calculated." ||
        echo "`date -R`: Dinucleotide abundace table cannot be calculated...") \
        >> {log} 2>&1

        (echo "`date -R`: Calculating nucleotide abundance table..." &&
        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {output.filt} \
        -k 1 \
        -o {output.nuc}  &&
        echo "`date -R`: Success! Nucleotide abundance table is calculated." ||
        echo "`date -R`: Nucleotide abundace table cannot be calculated...") \
        >> {log} 2>&1
        """