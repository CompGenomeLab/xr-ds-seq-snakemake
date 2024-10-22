
rule sort_filter:
    input:
        lambda w: input4rule(w, config["meta"], "sort_filter"),
    output:
        "results/{method}/{samples}/{samples}_{build}_sorted_chr.bed",
    params:
        filt=config["filter"],
    log:
        "logs/rule/sort_filter/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/sort_filter/{samples}_{build}_{method}.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Sorting and filtering bed file by chromosomes..." &&
        sort -k1,1 -k2,2n -k3,3n {input} |&
        egrep {params.filt} > {output} &&
        echo "`date -R`: Success! Bed file is filtered." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }} ) > {log} 2>&1
        """

rule length_distribution:
    input:
        rules.sort_filter.output
    output:
        "results/{method}/{samples}/{samples}_{build}_length_distribution.txt",
    log:
        "logs/rule/length_distribution/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/length_distribution/{samples}_{build}_{method}.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Calculating the read length distribution..." &&
        awk '{{print $3-$2}}' {input} |&
        sort -k1,1n |& 
        uniq -c |& 
        sed 's/\s\s*/ /g' |&
        awk '{{print $2"\\t"$1}}' > {output} &&
        echo "`date -R`: Success! Length distribution is calculated." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
        """

rule length_mode:
    input:
        bed=rules.sort_filter.output,
        ld=rules.length_distribution.output,
    output:
        "results/{method}/{samples}/{samples}_{build}_lengthMode.bed",
    log:
        "logs/rule/length_mode/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/length_mode/{samples}_{build}_{method}.benchmark.txt",
    shell:  
        """
        length="$(awk -v m=0 '{{if(m<$2){{m=$2;l=$1}}}}END{{print l}}' \
        {input.ld})" 

        (echo "`date -R`: Filtering the reads by the lengths..." &&
        awk -v num="$length" '{{ if ($3-$2 == num) {{ print }} }}' {input.bed} \
        > {output} &&
        echo "`date -R`: Success! Reads are filtered." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
        """

rule sep_strands:
    input:
       rules.sort_filter.output,
    output:
        plus="results/{method}/{samples}/{samples}_{build}_sorted_plus.bed",
        minus="results/{method}/{samples}/{samples}_{build}_sorted_minus.bed",
    log:
        "logs/rule/sep_strands/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/sep_strands/{samples}_{build}_{method}.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Separating plus stranded reads..." &&
        awk '{{if($6=="+"){{print}}}}' {input} > {output.plus} &&
        echo "`date -R`: Success! Reads are separated." || 
        {{ echo "`date -R`: Process failed..."; rm {output.plus}; exit 1; }}  )  > {log} 2>&1

        (echo "`date -R`: Separating minus stranded reads..." &&
        awk '{{if($6=="-"){{print}}}}' {input} > {output.minus} &&
        echo "`date -R`: Success! Reads are separated." || 
        {{ echo "`date -R`: Process failed..."; rm {output.minus}; exit 1; }}  )  >> {log} 2>&1
        """

rule sep_strands_sim:
    input:
       "results/{method}/{samples}/{samples}_{build}_{method}_sim.bed",
    output:
        plus="results/processed_files/{samples}_{build}_{method}_sim_plus.bed",
        minus="results/processed_files/{samples}_{build}_{method}_sim_minus.bed",
    log:
        "logs/rule/sep_strands_sim/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/sep_strands_sim/{samples}_{build}_{method}.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Separating plus stranded reads..." &&
        awk '{{if($6=="+"){{print}}}}' {input} > {output.plus} &&
        echo "`date -R`: Success! Reads are separated." || 
        {{ echo "`date -R`: Process failed..."; rm {output.plus}; exit 1; }}  )  > {log} 2>&1

        (echo "`date -R`: Separating minus stranded reads..." &&
        awk '{{if($6=="-"){{print}}}}' {input} > {output.minus} &&
        echo "`date -R`: Success! Reads are separated." || 
        {{ echo "`date -R`: Process failed..."; rm {output.minus}; exit 1; }}  )  >> {log} 2>&1
        """

rule reposition:
    input:
        plus=rules.sep_strands.output.plus,
        minus=rules.sep_strands.output.minus,
        index=rules.genome_indexing.output,
    output:
        plus=temp("results/{method}/{samples}/{samples}_{build}_sorted_plus_10.bed"),
        minus=temp("results/{method}/{samples}/{samples}_{build}_sorted_minus_10.bed"), 
    log:
        "logs/rule/reposition/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/reposition/{samples}_{build}_{method}.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Centering {input.plus} at damage site..." &&
        bedtools flank \
        -i {input.plus} \
        -g {input.index} \
        -l 6 \
        -r 0 |& 
        bedtools slop \
        -g {input.index} \
        -l 0 \
        -r 4 |& 
        awk '{{ if ($3-$2 == 10) {{ print }} }}' \
        > {output.plus} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.plus}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Centering {input.minus} at damage site..." &&
        bedtools flank \
        -i {input.minus} \
        -g {input.index} \
        -l 0 \
        -r 6 |& 
        bedtools slop \
        -g {input.index} \
        -l 4 \
        -r 0 |& 
        awk '{{ if ($3-$2 == 10) {{ print }} }}' \
        > {output.minus} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.minus}; exit 1; }}  ) >> {log} 2>&1
        """

rule bed2fasta_ds:
    input:
        plus=rules.reposition.output.plus,
        minus=rules.reposition.output.minus,
        genome=rules.genome_download.output,
    output:
        plus=temp("results/{method}/{samples}/{samples}_{build}_sorted_plus_10.fa"),
        minus=temp("results/{method}/{samples}/{samples}_{build}_sorted_minus_10.fa"),
        comb=temp("results/{method}/{samples}/{samples}_{build}_sorted_10.fa"),
        bed=temp("results/{method}/{samples}/{samples}_{build}_sorted_10.bed"),       
    log:
        "logs/rule/bed2fasta_ds/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/bed2fasta_ds/{samples}_{build}_{method}.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        (echo "`date -R`: Combine files..." &&
        cat {input.plus} {input.minus} > {output.bed} &&
        echo "`date -R`: Success! Files are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        
        (echo "`date -R`: Converting {output.bed} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {output.bed} \
        -fo {output.comb} \
        -s &&
        echo "`date -R`: Success! {output.bed} is converted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Converting {input.plus} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.plus} \
        -fo {output.plus} \
        -s &&
        echo "`date -R`: Success! {input.plus} is converted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Converting {input.minus} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.minus} \
        -fo {output.minus} \
        -s &&
        echo "`date -R`: Success! {input.minus} is converted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """

rule bed2fasta_xr:
    input:
        bed=rules.length_mode.output,
        genome=rules.genome_download.output,
    output:
        temp("results/{method}/{samples}/{samples}_{build}_lengthMode.fa"),
    log:
        "logs/rule/bed2fasta_xr/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/bed2fasta_xr/{samples}_{build}_{method}.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        (echo "`date -R`: Converting {input.bed} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.bed} \
        -fo {output} \
        -s &&
        echo "`date -R`: Success! {input.bed} is converted." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
        """

rule bed2fasta_input:
    input:
        bed=lambda w: input4rule(w, config["meta"], "bed2fasta_input"),
        genome=rules.genome_download.output,
    output:
        "results/input/{samples}/{samples}_{build}.fasta",
    log:
        "logs/rule/bed2fasta_input/{samples}_{build}.log",
    benchmark:
        "logs/rule/bed2fasta_input/{samples}_{build}.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        (echo "`date -R`: Converting {input.bed} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.bed} \
        -fo {output} \
        -s &&
        echo "`date -R`: Success! {input.bed} is converted." || 
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
        """

rule bed2fasta_ds_after_filt:
    input:
        plus="results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        minus="results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed",
        genome=rules.genome_download.output,
    output:
        comb=temp("results/{method}/{samples}/{samples}_{build}_sorted_filt_10.fa"),
        bed=temp("results/{method}/{samples}/{samples}_{build}_sorted_filt_10.bed"),       
    log:
        "logs/rule/bed2fasta_ds_after_filt/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/bed2fasta_ds_after_filt/{samples}_{build}_{method}.benchmark.txt",
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        (echo "`date -R`: Combine files..." &&
        cat {input.plus} {input.minus} > {output.bed} &&
        echo "`date -R`: Success! Files are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        
        (echo "`date -R`: Converting {output.bed} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {output.bed} \
        -fo {output.comb} \
        -s &&
        echo "`date -R`: Success! {output.bed} is converted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """

rule filtbyMotifs:
    input:
        "results/{method}/{samples}/{samples}_{build}_sorted_{strand}_10.fa",
    output:
        "results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_{strand}.bed",
    params:
        lambda w: getMotif(w.samples, config["meta"][w.samples]["product"]),
    log:
        "logs/rule/filtbyMotifs/{samples}_{build}_{method}_{strand}.log",
    benchmark:
        "logs/rule/filtbyMotifs/{samples}_{build}_{method}_{strand}.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Filtering by the given motif..." &&
        workflow/scripts/fa2bedByChoosingReadMotifs.py \
        -i {input} \
        -o {output} \
        -r {params} &&
        echo "`date -R`: Success! Filtering is done." ||
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) > {log} 2>&1
        """

rule comb_strands:
    input:
        plus="results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        minus="results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed",
    output:
        "results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines.bed",
    log:
        "logs/rule/comb_strands/{samples}_{build}_{method}.log",
    benchmark:
        "logs/rule/comb_strands/{samples}_{build}_{method}.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Combine files..." &&
        cat {input.plus} {input.minus} > {output} &&
        echo "`date -R`: Success! Files are combined." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1
        """


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
        sort -k1,1 -k2,2n -k3,3n {input.p_bed} | 
        bedtools genomecov \
        -i - \
        -g {input.ref_genome} \
        -bg \
        -scale $(echo {params.read} | awk '{{print 1000000/$1}}') \
        > {output.plus} &&
        echo "`date -R`: Success! Genome coverage is calculated." || 
        {{ echo "`date -R`: Process failed..."; rm {output.plus}; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Calculating genome coverage of {input.m_bed}..." &&
        sort -k1,1 -k2,2n -k3,3n {input.m_bed} | 
        bedtools genomecov \
        -i - \
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
        
rule copy_bed_ds:
    input:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed",
    output:
        plus="results/processed_files/{samples}_{build}_DS_plus.bed",
        minus="results/processed_files/{samples}_{build}_DS_minus.bed",
    log:
        "logs/rule/copy_bed_ds/{samples}_{build}.log",
    benchmark:
        "logs/rule/copy_bed_ds/{samples}_{build}.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Copy plus stranded reads..." &&
        cp {input.plus} {output.plus} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.plus}; exit 1; }}  )  > {log} 2>&1

        (echo "`date -R`: Copy minus stranded reads..." &&
        cp {input.minus} {output.minus} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.minus}; exit 1; }}  )  >> {log} 2>&1
        """

rule copy_bed_xr:
    input:
        plus="results/XR/{samples}/{samples}_{build}_sorted_plus.bed",
        minus="results/XR/{samples}/{samples}_{build}_sorted_minus.bed",
    output:
        plus="results/processed_files/{samples}_{build}_XR_plus.bed",
        minus="results/processed_files/{samples}_{build}_XR_minus.bed",
    log:
        "logs/rule/copy_bed_xr/{samples}_{build}.log",
    benchmark:
        "logs/rule/copy_bed_xr/{samples}_{build}.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Copy plus stranded reads..." &&
        cp {input.plus} {output.plus} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.plus}; exit 1; }}  )  > {log} 2>&1

        (echo "`date -R`: Copy minus stranded reads..." &&
        cp {input.minus} {output.minus} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; rm {output.minus}; exit 1; }}  )  >> {log} 2>&1
        """