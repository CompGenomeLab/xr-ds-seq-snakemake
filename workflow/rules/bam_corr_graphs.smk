rule bam_corr_graphs:
    input:
        npz="results/readCounts_{method}.npz",
    output:
        scatter=report("results/scatterplot_PearsonCorr_bigwigScores_{method}.png", 
                category="Correlation"),
        tab="results/PearsonCorr_bigwigScores_{method}.tab",
        heatmap=report("results/heatmap_SpearmanCorr_readCounts_{method}.png", 
                category="Correlation"),
        tab2="results/SpearmanCorr_readCounts_{method}.tab",
        pca=report("results/PCA_readCounts_{method}.png", 
                category="Correlation"),
    log:
        "logs/bam_corr_graphs_{method}.log",
    benchmark:
        "logs/bam_corr_graphs.benchmark_{method}.txt",
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
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1

        (echo "`date -R`: Plotting correlation (heatmap)..." &&
        plotCorrelation \
        -in {input.npz} \
        --corMethod spearman --skipZeros \
        --plotTitle "Spearman Correlation of Read Counts" \
        --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
        -o {output.heatmap} \
        --outFileCorMatrix {output.tab2} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1

        (echo "`date -R`: PCA analysis..." &&
        plotPCA -in {input.npz} \
        -o {output.pca} \
        -T "PCA of read counts" &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) >> {log} 2>&1
        """