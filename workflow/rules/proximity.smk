rule proximity_test:
    input:
        seeds="results/genes/seeds_{disease}.tsv",
        network="results/networks/string_bg_{disease}.tsv",
        expanded="results/genes/expanded_{disease}_rwr.tsv"
    output:
        significant="results/genes/significant_{disease}.tsv",
        meta="results/genes/significant_{disease}.tsv.meta.json"
    conda:
        "../../envs/proximity.yaml"
    script:
        "../scripts/proximity_test.py"
