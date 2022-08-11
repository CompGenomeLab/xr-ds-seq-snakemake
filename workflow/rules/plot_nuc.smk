
rule nucleotide_table:
    input:
        lambda w: input4nucTable(w, config["meta"]),
    output:
        dinuc=temp("results/{method}/{samples}/{samples}_{build}_sorted_dinucleotideTable.txt"),
        nuc=temp("results/{method}/{samples}/{samples}_{build}_sorted_nucleotideTable.txt"),
    log:
        "logs/rule/nucleotide_table/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/nucleotide_table/{samples}_{build}_{method}.benchmark.txt",
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

rule nucleotide_table_ds_sim:
    input:
        "results/{method}/{samples}/{samples}_{build}_ds_sim.fa",
    output:
        dinuc=temp("results/{method}/{samples}/{samples}_{build}_DS_sim_dinucleotideTable.txt"),
        nuc=temp("results/{method}/{samples}/{samples}_{build}_DS_sim_nucleotideTable.txt"),
    log:
        "logs/rule/nucleotide_table_ds_sim/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/nucleotide_table_ds_sim/{samples}_{build}_{method}.benchmark.txt",
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
        -o {output.nuc} &&
        echo "`date -R`: Success! Nucleotide abundance table is calculated." ||
        {{ echo "`date -R`: Nucleotide abundace table cannot be calculated..."; rm {output.nuc}; exit 1; }}  ) \
        >> {log} 2>&1
        """

rule nucleotide_table_xr_sim:
    input:
        "results/{method}/{samples}/{samples}_{build}_xr_sim.fa",
    output:
        dinuc=temp("results/{method}/{samples}/{samples}_{build}_XR_sim_dinucleotideTable.txt"),
        nuc=temp("results/{method}/{samples}/{samples}_{build}_XR_sim_nucleotideTable.txt"),
        filt=temp("results/{method}/{samples}/{samples}_{build}_XR_sim_filt.fa"),
    log:
        "logs/rule/nucleotide_table_xr_sim/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/nucleotide_table_xr_sim/{samples}_{build}_{method}.benchmark.txt",
    shell:  
        """
        len="$(awk 'BEGIN{{RS=">";ORS=""}}{{print length($2)"\\n"}}' {input} | sort | uniq -c | sort -k1,1n | tail -1 | awk '{{print $2}}')"

        ( echo "`date -R`: Filtering by most occurred read length..." &&
        awk -v num="$len" \
        'BEGIN{{RS=">";ORS=""}}length($2)==num{{print ">"$0}}' \
        {input} > {output.filt} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Filtering is failed..."; exit 1; }} ) \
        > {log} 2>&1
        
        ( echo "`date -R`: Calculating dinucleotide abundance table..." &&
        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {output.filt} \
        -k 2 \
        -o {output.dinuc} &&
        echo "`date -R`: Success! Dinucleotide abundance table is calculated." ||
        {{ echo "`date -R`: Dinucleotide abundace table cannot be calculated..."; 
        exit 1; }} ) >> {log} 2>&1

        ( echo "`date -R`: Calculating nucleotide abundance table..." &&
        workflow/scripts/fa2kmerAbundanceTable.py \
        -i {output.filt} \
        -k 1 \
        -o {output.nuc}  &&
        echo "`date -R`: Success! Nucleotide abundance table is calculated." ||
        {{ echo "`date -R`: Nucleotide abundace table cannot be calculated..."; \
        exit 1; }} ) >> {log} 2>&1
        """

rule plot_nuc:
    input:
        dinuc=rules.nucleotide_table.output.dinuc,
        nuc=rules.nucleotide_table.output.nuc,
    output:
        dinuc=report("results/processed_files/{samples}_{build}_{method}_dinucleotideTable.pdf", 
                    category="Nucleotide Content"),
        nuc=report("results/processed_files/{samples}_{build}_{method}_nucleotideTable.pdf", 
                    category="Nucleotide Content"),
    params:
        motif=lambda w: getDinuc(w.samples, config["meta"][w.samples]["product"]),
        name="{samples}_{build}_{method}",
    log:
        "logs/rule/plot_nuc/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/plot_nuc/{samples}_{build}_{method}.benchmark.txt",
    conda:
        "../envs/plot_nuc.yaml"
    shell:  
        """
        Rscript workflow/scripts/nucleotidePlot.r \
        -i {input.dinuc} \
        -k 2 \
        -s {params.name} \
        -f {params.motif} \
        -o {output.dinuc} \
        -l {log}

        Rscript workflow/scripts/nucleotidePlot.r \
        -i {input.nuc} \
        -k 1 \
        -s {params.name} \
        -o {output.nuc} \
        -l {log}
        """

rule plot_nuc_sim:
    input:
        dinuc="results/{method}/{samples}/{samples}_{build}_{method}_sim_dinucleotideTable.txt",
        nuc="results/{method}/{samples}/{samples}_{build}_{method}_sim_nucleotideTable.txt",
    output:
        dinuc=report("results/processed_files/{samples}_{build}_{method}_sim_dinucleotideTable.pdf", 
                    category="Nucleotide Content"),
        nuc=report("results/processed_files/{samples}_{build}_{method}_sim_nucleotideTable.pdf", 
                    category="Nucleotide Content"),
    params:
        motif=lambda w: getDinuc(w.samples, config["meta"][w.samples]["product"]),
        name="{samples}_{build}_{method}_sim",
    log:
        "logs/rule/plot_nuc_sim/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/plot_nuc_sim/{samples}_{build}_{method}.benchmark.txt",
    conda:
        "../envs/plot_nuc.yaml"
    shell:  
        """
        Rscript workflow/scripts/nucleotidePlot.r \
        -i {input.dinuc} \
        -k 2 \
        -s {params.name} \
        -f {params.motif} \
        -o {output.dinuc} \
        -l {log}

        Rscript workflow/scripts/nucleotidePlot.r \
        -i {input.nuc} \
        -k 1 \
        -s {params.name} \
        -o {output.nuc} \
        -l {log}
        """