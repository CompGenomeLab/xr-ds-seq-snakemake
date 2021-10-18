rule sra_se:
    output:
        "resources/samples/{samples}.fastq.gz", 
    params:
        srr=lambda w: getSRR(w.samples, config["srr"]["codes"], 
            config["sample"]),
        name="{samples}",
    log:
        "logs/{samples}/{samples}_se_sra.log",
    benchmark:
        "logs/{samples}/{samples}_se_sra.benchmark.txt",
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
            fasterq-dump \
            --threads {threads} \
            --progress $srr \
            -t resources/samples/ \
            -o resources/samples/${{srr}}.fastq &&
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
        srr=lambda w: getSRR(w.samples, config["srr"]["codes"], 
            config["sample"]),
        name="{samples}",
    log:
        "logs/{samples}/{samples}_pe_sra.log",
    benchmark:
        "logs/{samples}/{samples}_pe_sra.benchmark.txt",
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
            fasterq-dump \
            --threads {threads} \
            --progress $srr \
            -t resources/samples/ \
            -o resources/samples/${{srr}} &&
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