rule bam_correlation:
    input:
        lambda w: input4PCA(config["meta"], config["genome"]["build"]),
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