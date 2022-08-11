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