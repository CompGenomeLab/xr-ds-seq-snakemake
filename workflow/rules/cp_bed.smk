rule cp_bed_ds:
    input:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed",
    output:
        plus="results/processed_files/{samples}_{build}_DS_plus.bed",
        minus="results/processed_files/{samples}_{build}_DS_minus.bed",
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_cp_bed_ds.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_cp_bed_ds.benchmark.txt",
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

rule cp_bed_xr:
    input:
        plus="results/XR/{samples}/{samples}_{build}_sorted_plus.bed",
        minus="results/XR/{samples}/{samples}_{build}_sorted_minus.bed",
    output:
        plus="results/processed_files/{samples}_{build}_XR_plus.bed",
        minus="results/processed_files/{samples}_{build}_XR_minus.bed",
    log:
        "logs/rule/analysis/{samples}/{samples}_{build}_cp_bed_xr.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_{build}_cp_bed_xr.benchmark.txt",
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