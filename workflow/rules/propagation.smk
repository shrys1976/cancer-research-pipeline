rule run_rwr:
    input:
        seeds="results/genes/seeds_{disease}.tsv",
        network="results/networks/string_bg_{disease}.tsv"
    output:
        expanded="results/genes/expanded_{disease}_rwr.tsv",
        meta="results/genes/expanded_{disease}_rwr.tsv.meta.json"
    conda:
        "../../envs/propagation.yaml"
    script:
        "../scripts/rwr.py"
