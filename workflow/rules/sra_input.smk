
rule sra_se_input:
    output:
        "resources/input/{samples}.fastq.gz", 
    params:
        srr=lambda w: getSRR(w.samples, config["input"]["srr"]["codes"], 
            config["input"]["files"]),
        name="{samples}",
    log:
        "logs/rule/analysis/{samples}/{samples}_se_sra_input.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_se_sra_input.benchmark.txt",
    wildcard_constraints:
        samples='|'.join([x for x in config["input"]["files"]])
    conda:
        "../envs/sra.yaml"
    threads:
        6
    shell:
        """
        touch resources/input/{params.name}.fastq
        touch {log}

        srrList=$(echo {params.srr} | tr ":" "\\n")
        echo $srrList

        for srr in $srrList; do
            (echo "`date -R`: Downloading $srr files..." &&
            prefetch $srr \
            -O resources/samples/input/ &&
            vdb-validate resources/samples/input/$srr &&
            fastq-dump \
            resources/samples/input/$srr \
            --outdir resources/samples/input/ &&
            echo "`date -R`: Download is successful!" || 
            {{ echo "`date -R`: Process failed..."; exit 1; }} ) \
            >> {log} 2>&1 

            cat resources/input/${{srr}}.fastq >> resources/input/{params.name}.fastq

            rm resources/input/${{srr}}.fastq ; done

        (echo "`date -R`: Zipping srr file..." &&
        gzip resources/input/{params.name}.fastq &&
        echo "`date -R`: Zipping is successful!" || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }} ) \
        >> {log} 2>&1
        """

rule sra_pe_input:
    output:
        "resources/input/{samples}_1.fastq.gz", 
        "resources/input/{samples}_2.fastq.gz", 
    params:
        srr=lambda w: getSRR(w.samples, config["input"]["srr"]["codes"], 
            config["input"]["files"]),
        name="{samples}",
    log:
        "logs/rule/analysis/{samples}/{samples}_pe_sra_input.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_pe_sra_input.benchmark.txt",
    wildcard_constraints:
        samples='|'.join([x for x in config["input"]["files"]])
    conda:
        "../envs/sra.yaml"
    threads:
        6
    shell:
        """
        touch resources/samples/input/{params.name}_1.fastq
        touch resources/samples/input/{params.name}_2.fastq
        touch {log}
        echo "`date -R`: paired-end layout" >> {log}

        srrList=$(echo {params.srr} | tr ":" "\\n")
        echo $srrList >> {log}

        for srr in $srrList; do
            (echo "`date -R`: Downloading $srr files..." &&
            prefetch $srr \
            -O resources/samples/input/ &&
            vdb-validate resources/samples/input/$srr &&
            fastq-dump \
            resources/samples/input/$srr \
            --outdir resources/samples/input/ \
            --split-files &&
            echo "`date -R`: Download is successful!" || 
            {{ echo "`date -R`: Process failed..."; exit 1; }}  ) \
            > {log} 2>&1

            cat resources/samples/input/${{srr}}_1.fastq >> resources/samples/input/{params.name}_1.fastq
            cat resources/samples/input/${{srr}}_2.fastq >> resources/samples/input/{params.name}_2.fastq

            rm resources/samples/input/${{srr}}_1.fastq
            rm resources/samples/input/${{srr}}_2.fastq ; done

        (echo "`date -R`: Zipping srr file..." &&
        gzip resources/samples/input/{params.name}_1.fastq &&
        gzip resources/samples/input/{params.name}_2.fastq &&
        echo "`date -R`: Zipping is successful!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) \
        >> {log} 2>&1
        """