
rule simulation_ds_input:
    input:
        plus="results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        minus="results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed",
        genome=rules.genome_download.output,
        inpfile=lambda w: getInput(w.samples, config["meta"], config["genome"]["build"]),
    output:
        bed=temp("results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines.bed"),
        fa=temp("results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines.fa"),
        sim="results/{method}/{samples}/{samples}_{build}_ds_sim.fa",
        sam=temp("results/{method}/{samples}/{samples}_cutadapt_ds_sim_{build}.sam"),
        bam=temp("results/{method}/{samples}/{samples}_cutadapt_ds_sim_{build}.bam"),
        bam_sorted=temp("results/{method}/{samples}/{samples}_cutadapt_ds_sim_{build}_sorted.bam"),
        idx=temp("results/{method}/{samples}/{samples}_cutadapt_ds_sim_{build}_sorted.bam.bai"),
        simbed="results/{method}/{samples}/{samples}_{build}_ds_sim.bed",  
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_simulation_ds_input.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_simulation_ds_input.benchmark.txt",
    conda:
        "../envs/simulation.yaml"
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
        -fo {output.fa} \
        -s &&
        echo "`date -R`: Success! {output.bed} is converted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Simulating reads with input file..." &&
        boquila \
        --fasta {output.fa} \
        --inseqFasta \
        --inseq {input.inpfile} \
        --seed 1 \
        --sens 2 \
        > {output.sim} &&
        echo "`date -R`: Success! Simulation is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Aligning fasta file..." &&
        bowtie2 \
        -x {params.ref_genome} \
        -f {output.sim} -S {output.sam} &&
        echo "`date -R`: Success! Alignment is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Converting sam to bam..." &&
        samtools view -Sbh -o {output.bam} {output.sam} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Sorting (coordinates) bam file..." &&
        samtools sort {output.bam} > {output.bam_sorted} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Index bam file..." &&
        samtools index {output.bam_sorted} {output.idx} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Processing bam file..." && 
        samtools view -b {output.bam_sorted} |&
        bedtools bamtobed > {output.simbed} &&
        echo "`date -R`: Success! Bam file converted to bed format." || 
        {{ echo "`date -R`: Process failed..."; rm {output.simbed}; exit 1; }}  ) >> {log} 2>&1
        """

rule simulation_xr_input:
    input:
        bed=rules.sort_filter.output,
        genome=rules.genome_download.output,
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2", 
        inpfile=lambda w: getInput(w.samples, config["meta"], config["genome"]["build"]),
    output:
        fa=temp("results/{method}/{samples}/{samples}_{build}_sorted_chr.fa"),
        sim="results/{method}/{samples}/{samples}_{build}_xr_sim.fa",
        sam=temp("results/{method}/{samples}/{samples}_cutadapt_xr_sim_{build}.sam"),
        bam=temp("results/{method}/{samples}/{samples}_cutadapt_xr_sim_{build}.bam"),
        bam_sorted=temp("results/{method}/{samples}/{samples}_cutadapt_xr_sim_{build}_sorted.bam"),
        idx=temp("results/{method}/{samples}/{samples}_cutadapt_xr_sim_{build}_sorted.bam.bai"),
        simbed="results/{method}/{samples}/{samples}_{build}_xr_sim.bed",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_simulation_xr_input.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_simulation_xr_input.benchmark.txt",
    conda:
        "../envs/simulation.yaml"
    shell:
        """
        (echo "`date -R`: Converting {input.bed} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.bed} \
        -fo {output.fa} \
        -s &&
        echo "`date -R`: Success! {input.bed} is converted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Simulating reads with input file..." &&
        boquila \
        --fasta {output.fa} \
        --inseqFasta \
        --inseq {input.inpfile} \
        --seed 1 \
        --sens 2 \
        > {output.sim} &&
        echo "`date -R`: Success! Simulation is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Aligning fasta file..." &&
        bowtie2 \
        -x {params.ref_genome} \
        -f {output.sim} -S {output.sam} &&
        echo "`date -R`: Success! Alignment is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Converting sam to bam..." &&
        samtools view -Sbh -o {output.bam} {output.sam} &&
        echo "`date -R`: Success! Conversion is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Sorting (coordinates) bam file..." &&
        samtools sort {output.bam} > {output.bam_sorted} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Index bam file..." &&
        samtools index {output.bam_sorted} {output.idx} &&
        echo "`date -R`: Success! Bam file is sorted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Processing bam file..." && 
        samtools view -b {output.bam_sorted} |&
        bedtools bamtobed > {output.simbed} &&
        echo "`date -R`: Success! Bam file converted to bed format." || 
        {{ echo "`date -R`: Process failed..."; rm {output.simbed}; exit 1; }}  ) >> {log} 2>&1
        """

'''
rule simulation_ds:
    input:
        plus="results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        minus="results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed",
        genome="resources/ref_genomes/{build}/genome_{build}.fa",
        inpfile=lambda w: getInput(w.samples, config["meta"], config["genome"]["build"]),
    output:
        bed=temp("results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines.bed"),
        fa=temp("results/{method}/{samples}/{samples}_{build}_sorted_ds_dipyrimidines.fa"),
        sim="results/{method}/{samples}/{samples}_{build}_ds_sim.fa",
        simbed="results/{method}/{samples}/{samples}_{build}_ds_sim.bed",  
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_simulation_ds.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_simulation_ds.benchmark.txt",
    conda:
        "../envs/simulation.yaml"
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
        -fo {output.fa} \
        -s &&
        echo "`date -R`: Success! {output.bed} is converted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Simulating reads according to reference genome..." &&
        boquila \
        --fasta {output.fa} \
        --bed {output.simbed} \
        --ref {input.genome} \
        --seed 1 \
        --sens 2 \
        --regions {input.inpfile} \
        > {output.sim} &&
        echo "`date -R`: Success! Simulation is done." || 
        {{ echo "`date -R`: Process failed..."; rm {output.simbed}; exit 1; }}  ) >> {log} 2>&1
        """

rule simulation_xr:
    input:
        bed="results/{method}/{samples}/{samples}_{build}_sorted_chr.bed",
        genome="resources/ref_genomes/{build}/genome_{build}.fa",
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2", 
        inpfile=lambda w: getInput(w.samples, config["meta"], config["genome"]["build"]),
    output:
        fa=temp("results/{method}/{samples}/{samples}_{build}_sorted_chr.fa"),
        sim="results/{method}/{samples}/{samples}_{build}_xr_sim.fa",
        simbed="results/{method}/{samples}/{samples}_{build}_xr_sim.bed",
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_simulation_xr.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_{method}_simulation_xr.benchmark.txt",
    conda:
        "../envs/simulation.yaml"
    shell:
        """
        (echo "`date -R`: Converting {input.bed} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.bed} \
        -fo {output.fa} \
        -s &&
        echo "`date -R`: Success! {input.bed} is converted." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Simulating reads according to reference genome..." &&
        boquila \
        --fasta {output.fa} \
        --bed {output.simbed} \
        --ref {input.genome} \
        --seed 1 \
        --sens 2 \
        --regions {input.inpfile} \
        > {output.sim} &&
        echo "`date -R`: Success! Simulation is done." || 
        {{ echo "`date -R`: Process failed..."; rm {output.simbed}; exit 1; }}  ) >> {log} 2>&1
        """
'''