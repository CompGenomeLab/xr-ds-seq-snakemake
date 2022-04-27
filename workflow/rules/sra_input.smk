
rule sra_se_input:
    output:
        "resources/input/{samples}.fastq.gz", 
    params:
        name="{samples}",
    log:
        "logs/rule/analysis/{samples}/{samples}_se_sra_input.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_se_sra_input.benchmark.txt",
    conda:
        "../envs/sra.yaml"
    threads:
        6
    shell:
        """
        touch resources/input/{params.name}.fastq
        touch {log}

        srrList=$(echo {params.name} | tr ":" "\\n")
        echo "`date -R`: $srrList" >> {log}

        for srr in $srrList; do
            (echo "`date -R`: Downloading $srr files..." &&
            prefetch $srr \
            -O resources/input/ &&
            vdb-validate resources/input/$srr &&
            fastq-dump \
            resources/input/${{srr}}/${{srr}}.sra \
            --outdir resources/input/ &&
            echo "`date -R`: Download is successful!" || 
            {{ echo "`date -R`: Process failed..."; exit 1; }} ) \
            >> {log} 2>&1 

            cat resources/input/${{srr}}.fastq >> resources/input/{params.name}.fastq ; done

        (echo "`date -R`: Zipping srr file..." &&
        gzip resources/input/{params.name}.fastq &&
        echo "`date -R`: Zipping is successful!" || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }} ) \
        >> {log} 2>&1
        """

rule sra_pe_input:
    output:
        read1="resources/input/{samples}_1.fastq.gz", 
        read2="resources/input/{samples}_2.fastq.gz", 
    params:
        name="{samples}",
    log:
        "logs/rule/analysis/{samples}/{samples}_pe_sra_input.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_pe_sra_input.benchmark.txt",
    conda:
        "../envs/sra.yaml"
    threads:
        6
    shell:
        """
        touch resources/input/{params.name}_1.fastq
        touch resources/input/{params.name}_2.fastq
        touch {log}
        echo "`date -R`: paired-end layout" >> {log}

        srrList=$(echo {params.name} | tr ":" "\\n")
        echo "`date -R`: $srrList" >> {log}

        for srr in $srrList; do
            (echo "`date -R`: Downloading $srr files..." &&
            prefetch $srr \
            -O resources/input/ &&
            vdb-validate resources/input/$srr &&
            fastq-dump \
            resources/input/${{srr}}/${{srr}}.sra \
            --outdir resources/input/ \
            --split-files &&
            echo "`date -R`: Download is successful!" || 
            {{ echo "`date -R`: Process failed..."; exit 1; }}  ) \
            >> {log} 2>&1

            cat resources/input/${{srr}}_1.fastq >> resources/input/{params.name}_1.fastq
            cat resources/input/${{srr}}_2.fastq >> resources/input/{params.name}_2.fastq ; done

        (echo "`date -R`: Zipping srr file..." &&
        gzip resources/input/{params.name}_1.fastq &&
        gzip resources/input/{params.name}_2.fastq &&
        echo "`date -R`: Zipping is successful!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) \
        >> {log} 2>&1
        """