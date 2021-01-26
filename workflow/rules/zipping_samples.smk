
rule zip_samples:
    input:
        "resources/samples/{samples}"
    output:
        "resources/samples/{samples}.gz"
    shell:
        "gzip {input}"