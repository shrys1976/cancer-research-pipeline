rule run_enrichment:
    input:
        string_leiden="results/modules/modules_{disease}_leiden.tsv",
        string_louvain="results/modules/modules_{disease}_louvain.tsv",
        disease="results/modules/modules_diseasePPI_{disease}.tsv"
    output:
        enrichment="results/enrichment/enrichment_{disease}.tsv",
        meta="results/enrichment/enrichment_{disease}.tsv.meta.json"
    conda:
        "../../envs/enrichment.yaml"
    script:
        "../scripts/enrichment.py"
