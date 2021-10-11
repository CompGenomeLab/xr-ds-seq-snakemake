
rule sra_se_input:
    output:
        "resources/input/{samples}.fastq.gz", 
    params:
        srr=lambda w: getSRR(w.samples, config["input"]["srr"]["codes"], 
            config["input"]["files"]),
        name="{samples}",
    log:
        "logs/{samples}/{samples}_se_sra_input.log",
    benchmark:
        "logs/{samples}/{samples}_se_sra_input.benchmark.txt",
    wildcard_constraints:
        samples='|'.join([x for x in config["input"]["files"]])
    conda:
        "../envs/sra.yaml"
    threads:
        8
    shell:
        """
        touch resources/input/{params.name}.fastq
        touch {log}

        srrList=$(echo {params.srr} | tr ":" "\\n")
        echo $srrList

        for srr in $srrList; do
            (echo "`date -R`: Downloading $srr files..." &&
            fasterq-dump \
            --threads {threads} \
            --progress $srr \
            -t resources/input/ \
            -o resources/input/${{srr}}.fastq &&
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
        "logs/{samples}/{samples}_pe_sra_input.log",
    benchmark:
        "logs/{samples}/{samples}_pe_sra_input.benchmark.txt",
    wildcard_constraints:
        samples='|'.join([x for x in config["input"]["files"]])
    conda:
        "../envs/sra.yaml"
    threads:
        8
    shell:
        """
        touch resources/input/{params.name}_1.fastq
        touch resources/input/{params.name}_2.fastq
        touch {log}

        srrList=$(echo {params.srr} | tr ":" "\\n")
        echo $srrList

        for srr in $srrList; do

            (echo "`date -R`: Downloading SRR files..." &&
            fasterq-dump \
            --threads {threads} \
            --progress {params.srr} \
            -t resources/input/ \
            -o resources/input/${{srr}} &&
            echo "`date -R`: Download is successful!" || 
            {{ echo "`date -R`: Process failed..."; exit 1; }} ) \
            >> {log} 2>&1

            cat resources/input/${{srr}}_1.fastq >> resources/input/{params.name}_1.fastq
            cat resources/input/${{srr}}_2.fastq >> resources/input/{params.name}_2.fastq

            rm resources/input/${{srr}}_1.fastq 
            rm resources/input/${{srr}}_2.fastq ; done

        (echo "`date -R`: Zipping srr file..." &&
        gzip resources/input/{params.name}*.fastq &&
        echo "`date -R`: Zipping is successful!" || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }} ) \
        >> {log} 2>&1
        """