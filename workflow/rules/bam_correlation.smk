rule bam_correlation:
    input:
        lambda w: input4PCA(config["meta"], config["genome"]["build"]),
    output:
        out="results/readCounts.npz",
        raw_out="results/readCounts.tab",
    log:
        "logs/rule/analysis/bam_correlation.log",
    benchmark:
        "logs/rule/analysis/bam_correlation.benchmark.txt",
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
