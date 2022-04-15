rule sra_se:
    output:
        "resources/samples/{samples}.fastq.gz", 
    params:
        srr=lambda w: config["meta"][w.samples]["srr_id"],
        name="{samples}",
    log:
        "logs/rule/analysis/{samples}/{samples}_se_sra.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_se_sra.benchmark.txt",
    wildcard_constraints:
        samples='|'.join([x for x in config["sample"]])
    conda:
        "../envs/sra.yaml"
    threads:
        6
    shell:
        """
        touch resources/samples/{params.name}.fastq
        touch {log}

        srrList=$(echo {params.srr} | tr ":" "\\n")
        echo $srrList

        for srr in $srrList; do
            (echo "`date -R`: Downloading $srr files..." &&
            prefetch $srr \
            -O resources/samples/ &&
            vdb-validate resources/samples/$srr &&
            fastq-dump \
            resources/samples/$srr \
            --outdir resources/samples/ \
            echo "`date -R`: Download is successful!" || 
            {{ echo "`date -R`: Process failed..."; exit 1; }} ) \
            >> {log} 2>&1 

            cat resources/samples/${{srr}}.fastq >> resources/samples/{params.name}.fastq

            rm resources/samples/${{srr}}.fastq ; done

        (echo "`date -R`: Zipping srr file..." &&
        gzip resources/samples/{params.name}.fastq &&
        echo "`date -R`: Zipping is successful!" || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }} ) \
        >> {log} 2>&1
        """

rule sra_pe:
    output:
        "resources/samples/{samples}_1.fastq.gz", 
        "resources/samples/{samples}_2.fastq.gz", 
    params:
        srr=lambda w: config["meta"][w.samples]["srr_id"],
        name="{samples}",
    log:
        "logs/rule/analysis/{samples}/{samples}_pe_sra.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_pe_sra.benchmark.txt",
    wildcard_constraints:
        samples='|'.join([x for x in config["sample"]])
    conda:
        "../envs/sra.yaml"
    threads:
        6
    shell:
        """
        touch resources/samples/{params.name}_1.fastq
        touch resources/samples/{params.name}_2.fastq
        touch {log}
        echo "`date -R`: paired-end layout" >> {log}

        srrList=$(echo {params.srr} | tr ":" "\\n")
        echo $srrList >> {log}

        for srr in $srrList; do
            (echo "`date -R`: Downloading $srr files..." &&
            prefetch $srr \
            -O resources/samples/ &&
            vdb-validate resources/samples/$srr &&
            fastq-dump \
            resources/samples/$srr \
            --outdir resources/samples/ \
            --split-files &&
            echo "`date -R`: Download is successful!" || 
            {{ echo "`date -R`: Process failed..."; exit 1; }}  ) \
            > {log} 2>&1

            cat resources/samples/${{srr}}_1.fastq >> resources/samples/{params.name}_1.fastq
            cat resources/samples/${{srr}}_2.fastq >> resources/samples/{params.name}_2.fastq

            rm resources/samples/${{srr}}_1.fastq
            rm resources/samples/${{srr}}_2.fastq ; done

        (echo "`date -R`: Zipping srr file..." &&
        gzip resources/samples/{params.name}_1.fastq &&
        gzip resources/samples/{params.name}_2.fastq &&
        echo "`date -R`: Zipping is successful!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) \
        >> {log} 2>&1
        """