
rule plot_length:
    input:
        rules.length_distribution.output
    output:
        report("results/processed_files/{samples}_{build}_{method}_length_distribution.pdf", 
                category="Length Distribution"),
    params:
        "{samples}_{build}_{method}",
    log:
        "logs/rule/plot_length/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/plot_length/{samples}_{build}_{method}.benchmark.txt",
    conda:
        "../envs/plot_nuc.yaml"
    shell:  
        """
        Rscript workflow/scripts/lengthDistPlot.r \
        -i {input} \
        -s {params} \
        -o {output} \
        -l {log}
        """

rule nucleotide_table:
    input:
        lambda w: input4rule(w, config["meta"], "nucleotide_table"),
    output:
        dinuc="results/{method}/{samples}/{samples}_{build}_sorted_dinucleotideTable.txt",
        nuc="results/{method}/{samples}/{samples}_{build}_sorted_nucleotideTable.txt",
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

rule nucleotide_table_after_filt:
    input:
        lambda w: input4rule(w, config["meta"], "nucleotide_table", filtered = True),
    output:
        dinuc="results/{method}/{samples}/{samples}_{build}_sorted_filt_dinucleotideTable.txt",
        nuc="results/{method}/{samples}/{samples}_{build}_sorted_filt_nucleotideTable.txt",
    log:
        "logs/rule/nucleotide_table_after_filt/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/nucleotide_table_after_filt/{samples}_{build}_{method}.benchmark.txt",
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
        dinuc="results/{method}/{samples}/{samples}_{build}_DS_sim_dinucleotideTable.txt",
        nuc="results/{method}/{samples}/{samples}_{build}_DS_sim_nucleotideTable.txt",
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
        dinuc="results/{method}/{samples}/{samples}_{build}_XR_sim_dinucleotideTable.txt",
        nuc="results/{method}/{samples}/{samples}_{build}_XR_sim_nucleotideTable.txt",
        filt="results/{method}/{samples}/{samples}_{build}_XR_sim_filt.fa",
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

rule plot_nuc_after_filt:
    input:
        dinuc=rules.nucleotide_table_after_filt.output.dinuc,
        nuc=rules.nucleotide_table_after_filt.output.nuc,
    output:
        dinuc=report("results/processed_files/{samples}_{build}_{method}_filt_dinucleotideTable.pdf", 
                    category="Nucleotide Content"),
        nuc=report("results/processed_files/{samples}_{build}_{method}_filt_nucleotideTable.pdf", 
                    category="Nucleotide Content"),
    params:
        motif=lambda w: getDinuc(w.samples, config["meta"][w.samples]["product"]),
        name="{samples}_{build}_{method}_filt",
    log:
        "logs/rule/plot_nuc_after_filt/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/plot_nuc_after_filt/{samples}_{build}_{method}.benchmark.txt",
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

rule bam_correlation:
    input:
        lambda w: input4rule(w, config["meta"], "bam_correlation"),
    output:
        out="results/readCounts.npz",
        raw_out="results/readCounts.tab",
    log:
        "logs/rule/bam_correlation/bam_correlation.log",
    benchmark:
        "logs/rule/bam_correlation/bam_correlation.benchmark.txt",
    conda:
        "../envs/bam_correlation.yaml"
    shell:  
        """
        (echo "`date -R`: MultiBam summary..." &&
        multiBamSummary bins \
        --bamfiles {input} \
        --minMappingQuality 20 \
        -out {output.out} --outRawCounts {output.raw_out} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.out}; exit 1; }} ) > {log} 2>&1
        """

rule plot_bam_corr:
    input:
        npz=rules.bam_correlation.output.out
    output:
        scatter=report("results/processed_files/scatterplot_PearsonCorr_bigwigScores.pdf", 
                category="Correlation"),
        tab="results/processed_files/PearsonCorr_bigwigScores.tab",
        heatmap=report("results/processed_files/heatmap_SpearmanCorr_readCounts.pdf", 
                category="Correlation"),
        tab2="results/processed_files/SpearmanCorr_readCounts.tab",
        pca=report("results/processed_files/PCA_readCounts.pdf", 
                category="Correlation"),
    log:
        "logs/rule/plot_bam_corr/plot_bam_corr.log",
    benchmark:
        "logs/rule/plot_bam_corr/plot_bam_corr.benchmark.txt",
    conda:
        "../envs/bam_correlation.yaml"
    shell:  
        """
        (echo "`date -R`: Plotting correlation (scatter)..." &&
        plotCorrelation \
        -in {input.npz} \
        --corMethod pearson --skipZeros \
        --plotTitle "Pearson Correlation of Bam Files" \
        --whatToPlot scatterplot \
        -o {output.scatter} \
        --outFileCorMatrix {output.tab} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.scatter}; exit 1; }} ) > {log} 2>&1

        (echo "`date -R`: Plotting correlation (heatmap)..." &&
        plotCorrelation \
        -in {input.npz} \
        --corMethod spearman --skipZeros \
        --plotTitle "Spearman Correlation of Read Counts" \
        --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
        -o {output.heatmap} \
        --outFileCorMatrix {output.tab2} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.heatmap}; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: PCA analysis..." &&
        plotPCA -in {input.npz} \
        -o {output.pca} \
        -T "PCA of read counts" &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.pca}; exit 1; }} ) >> {log} 2>&1
        """