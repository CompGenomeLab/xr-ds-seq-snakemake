
rule length_mode:
    input:
        bed="results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_chr.bed",
        ld="results/{dir}/{samples}{v}/{samples}_cutadapt_length_distribution.txt",
    output:
        "results/{dir}/{samples}{v}/{samples}_cutadapt_sorted_lengthMode.bed",
    log:
        "results/{dir}/{samples}{v}/log/length_mode.log"
    benchmark:
        "results/{dir}/{samples}{v}/log/length_mode.benchmark.txt",
    shell:  
        """
        length="$(awk -v m=0 '{{if(m<$2){{m=$2;l=$1}}}}END{{print l}}' \
        {input.ld})"

        awk -v num="$length" '{{ if ($3-$2 == num) {{ print }} }}' {input.bed} \
        > {output}
        """


