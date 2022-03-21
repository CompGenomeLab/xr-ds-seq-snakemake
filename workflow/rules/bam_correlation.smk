rule bam_correlation:
    input:
        lambda w: input4PCA(config["sample"], config["srr"]["enabled"], 
            config["srr"]["codes"], config["build"], config["method"]),
    output:
        out="results/readCounts_{method}.npz",
        raw_out="results/readCounts_{method}.tab",
    log:
        "logs/rule/analysis/bam_correlation_{method}.log",
    benchmark:
        "logs/rule/analysis/bam_correlation_{method}.benchmark.txt",
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
