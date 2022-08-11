
rule genomecov_ds:
    input:
        p_bed="results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        m_bed="results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed",
        ref_genome=rules.genome_indexing.output,
    output:
        plus=temp("results/{method}/{samples}/{samples}_{build}_DS_sorted_plus.bdg"),
        minus=temp("results/{method}/{samples}/{samples}_{build}_DS_sorted_minus.bdg"),
    params:
        read=lambda w, input: mappedReads(input[0], input[1]),
    log:
        "logs/rule/genomecov_ds/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/genomecov_ds/{samples}_{build}_{method}.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Calculating genome coverage of {input.p_bed}..." &&
        bedtools genomecov \
        -i {input.p_bed} \
        -g {input.ref_genome} \
        -bg \
        -scale $(echo {params.read} | awk '{{print 1000000/$1}}') \
        > {output.plus} &&
        echo "`date -R`: Success! Genome coverage is calculated." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Calculating genome coverage of {input.m_bed}..." &&
        bedtools genomecov \
        -i {input.m_bed} \
        -g {input.ref_genome} \
        -bg \
        -scale $(echo {params.read} | awk '{{print 1000000/$1}}') \
        > {output.minus} &&
        echo "`date -R`: Success! Genome coverage is calculated." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1
        """

rule genomecov_xr:
    input:
        p_bed="results/{method}/{samples}/{samples}_{build}_sorted_plus.bed",
        m_bed="results/{method}/{samples}/{samples}_{build}_sorted_minus.bed",
        ref_genome=rules.genome_indexing.output,
    output:
        plus=temp("results/{method}/{samples}/{samples}_{build}_XR_sorted_plus.bdg"),
        minus=temp("results/{method}/{samples}/{samples}_{build}_XR_sorted_minus.bdg"),
    params:
        read=lambda w, input: mappedReads(input[0], input[1]),
    log:
        "logs/rule/genomecov_xr/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/genomecov_xr/{samples}_{build}_{method}.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Calculating genome coverage of {input.p_bed}..." &&
        bedtools genomecov \
        -i {input.p_bed} \
        -g {input.ref_genome} \
        -bg \
        -scale $(echo {params.read} | awk '{{print 1000000/$1}}') \
        > {output.plus} &&
        echo "`date -R`: Success! Genome coverage is calculated." || 
        {{ echo "`date -R`: Process failed..."; rm {output.plus}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Calculating genome coverage of {input.m_bed}..." &&
        bedtools genomecov \
        -i {input.m_bed} \
        -g {input.ref_genome} \
        -bg \
        -scale $(echo {params.read} | awk '{{print 1000000/$1}}') \
        > {output.minus} &&
        echo "`date -R`: Success! Genome coverage is calculated." || 
        {{ echo "`date -R`: Process failed..."; rm {output.minus}; exit 1; }}  ) >> {log} 2>&1
        """

rule genomecov_sim:
    input:
        p_bed="results/processed_files/{samples}_{build}_{method}_sim_plus.bed",
        m_bed="results/processed_files/{samples}_{build}_{method}_sim_minus.bed",
        ref_genome=rules.genome_indexing.output,
    output:
        plus=temp("results/{method}/{samples}/{samples}_{build}_{method}_sim_plus.bdg"),
        minus=temp("results/{method}/{samples}/{samples}_{build}_{method}_sim_minus.bdg"),
    params:
        read=lambda w, input: mappedReads(input[0], input[1]),
    log:
        "logs/rule/genomecov_sim/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/genomecov_sim/{samples}_{build}_{method}.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Calculating genome coverage of {input.p_bed}..." &&
        bedtools genomecov \
        -i {input.p_bed} \
        -g {input.ref_genome} \
        -bg \
        -scale $(echo {params.read} | awk '{{print 1000000/$1}}') \
        > {output.plus} &&
        echo "`date -R`: Success! Genome coverage is calculated." || 
        {{ echo "`date -R`: Process failed..."; rm {output.plus}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Calculating genome coverage of {input.m_bed}..." &&
        bedtools genomecov \
        -i {input.m_bed} \
        -g {input.ref_genome} \
        -bg \
        -scale $(echo {params.read} | awk '{{print 1000000/$1}}') \
        > {output.minus} &&
        echo "`date -R`: Success! Genome coverage is calculated." || 
        {{ echo "`date -R`: Process failed..."; rm {output.minus}; exit 1; }}  ) >> {log} 2>&1
        """
        
rule bedGraphToBigWig:
    input:
        bdg="results/{method}/{samples}/{samples}_{build}_{method}_sorted_{strand}.bdg",
        index=rules.genome_indexing.output,
    output:
        sort=temp("results/{method}/{samples}/{samples}_{build}_{method}_resorted_{strand}.bdg"),
        bw="results/processed_files/{samples}_{build}_{method}_{strand}.bw",
    log:
        "logs/rule/bedGraphToBigWig/{samples}_{build}_{method}_{strand}.log",
    benchmark:
        "logs/rule/bedGraphToBigWig/{samples}_{build}_{method}_{strand}.benchmark.txt",
    conda:
        "../envs/bedGraphToBigWig.yaml"
    shell:  
        """
        (echo "`date -R`: Converting bedGraph to bigWig..." &&
        LC_COLLATE=C sort -k1,1 -k2,2n {input.bdg} > {output.sort} &&
        bedGraphToBigWig {output.sort} {input.index} {output.bw} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
        """
    
rule bedGraphToBigWig_sim:
    input:
        bdg="results/{method}/{samples}/{samples}_{build}_{method}_sim_{strand}.bdg",
        index=rules.genome_indexing.output,
    output:
        sort=temp("results/{method}/{samples}/{samples}_{build}_{method}_sim_resorted_{strand}.bdg"),
        bw="results/processed_files/{samples}_{build}_{method}_sim_{strand}.bw", 
    log:
        "logs/rule/bedGraphToBigWig_sim/{samples}_{build}_{method}_{strand}.log",
    benchmark:
        "logs/rule/bedGraphToBigWig_sim/{samples}_{build}_{method}_{strand}.benchmark.txt",
    conda:
        "../envs/bedGraphToBigWig.yaml"
    shell:  
        """
        (echo "`date -R`: Converting bedGraph to bigWig..." &&
        LC_COLLATE=C sort -k1,1 -k2,2n {input.bdg} > {output.sort} &&
        bedGraphToBigWig {output.sort} {input.index} {output.bw} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
        """