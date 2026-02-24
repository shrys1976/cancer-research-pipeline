rule fetch_string_background:
    input:
        seeds="results/genes/seeds_{disease}.tsv"
    output:
        network="results/networks/string_bg_{disease}.tsv",
        meta="results/networks/string_bg_{disease}.tsv.meta.json"
    conda:
        "../../envs/string.yaml"
    script:
        "../scripts/fetch_string.py"
