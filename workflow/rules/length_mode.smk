
rule length_mode:
    input:
        bed="results/{samples}/{samples}_{build}_sorted_chr.bed",
        ld="results/{samples}/{samples}_{build}_length_distribution.txt",
    output:
        temp("results/{samples}/{samples}_{build}_lengthMode.bed"),
    log:
        "logs/{samples}/{samples}_{build}_length_mode.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_length_mode.benchmark.txt",
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