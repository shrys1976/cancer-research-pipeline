rule query_psicquic:
    output:
        network="results/networks/disease_{disease}.tsv",
        meta="results/networks/disease_{disease}.tsv.meta.json"
    conda:
        "../../envs/psicquic.yaml"
    script:
        "../scripts/psicquic_query.py"
