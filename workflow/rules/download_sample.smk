rule sra_se:
    output:
        "resources/samples/{samples}.fastq.gz", 
    params:
        srr=lambda w: config["meta"][w.samples]["srr_id"],
        name="{samples}",
    log:
        "logs/rule/sra_se/{samples}.log",
    benchmark:
        "logs/rule/sra_se/{samples}.benchmark.txt",
    wildcard_constraints:
        samples='|'.join([x for x in config["meta"].keys()])
    conda:
        "../envs/sra.yaml"
    threads:
        6
    shell:
        """
        touch resources/samples/{params.name}.fastq
        touch {log}

        srrList=$(echo {params.srr} | tr ":" "\\n")
        echo $srrList >> {log}

        for srr in $srrList; do
            (echo "`date -R`: Downloading $srr files..." &&
            prefetch $srr \
            -O resources/samples/ &&
            vdb-validate resources/samples/$srr &&
            fastq-dump \
            resources/samples/${{srr}}/${{srr}}.sra \
            --outdir resources/samples/ &&
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
        read1="resources/samples/{samples}_1.fastq.gz", 
        read2="resources/samples/{samples}_2.fastq.gz", 
    params:
        srr=lambda w: config["meta"][w.samples]["srr_id"],
        name="{samples}",
    log:
        "logs/rule/sra_pe/{samples}.log",
    benchmark:
        "logs/rule/sra_pe/{samples}.benchmark.txt",
    wildcard_constraints:
        samples='|'.join([x for x in config["meta"].keys()])
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
            resources/samples/${{srr}}/${{srr}}.sra \
            --outdir resources/samples/ \
            --split-files &&
            echo "`date -R`: Download is successful!" || 
            {{ echo "`date -R`: Process failed..."; exit 1; }}  ) \
            >> {log} 2>&1

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

rule sra_se_input:
    output:
        "resources/input/{samples}.fastq.gz", 
    params:
        name="{samples}",
    log:
        "logs/rule/sra_se_input/{samples}.log",
    benchmark:
        "logs/rule/sra_se_input/{samples}.benchmark.txt",
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
        "logs/rule/sra_pe_input/{samples}.log",
    benchmark:
        "logs/rule/sra_pe_input/{samples}.benchmark.txt",
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

rule rename_raw:
    input:
        r1=lambda w: getPaired(w.samples, "forward", "resources/samples/"),
        r2=lambda w: getPaired(w.samples, "reverse", "resources/samples/"),
    output:
        r1="resources/samples/{samples}_1.fastq.gz", 
        r2="resources/samples/{samples}_2.fastq.gz",
    log:
        "logs/rule/rename_raw/{samples}.log",
    benchmark:
        "logs/rule/rename_raw/{samples}.benchmark.txt",
    shell:
        """
        mv {input.r1} {output.r1}

        mv {input.r2} {output.r2}
        """

rule rename_raw_input:
    input:
        r1=lambda w: getPaired(w.samples, "forward", "resources/input/"),
        r2=lambda w: getPaired(w.samples, "reverse", "resources/input/"),
    output:
        r1="resources/input/{samples}_1.fastq.gz", 
        r2="resources/input/{samples}_2.fastq.gz",
    log:
        "logs/rule/rename_raw_input/{samples}.log",
    benchmark:
        "logs/rule/rename_raw_input/{samples}.benchmark.txt",
    shell:
        """
        mv {input.r1} {output.r1}

        mv {input.r2} {output.r2}
        """