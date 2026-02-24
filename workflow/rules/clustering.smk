rule cluster_string_modules:
    input:
        network="results/networks/string_bg_{disease}.tsv",
        significant="results/genes/significant_{disease}.tsv"
    output:
        modules="results/modules/modules_{disease}_{method}.tsv",
        meta="results/modules/modules_{disease}_{method}.tsv.meta.json"
    params:
        mode="string_significant"
    wildcard_constraints:
        method="leiden|louvain"
    conda:
        "../../envs/clustering.yaml"
    script:
        "../scripts/cluster_graph.py"


rule cluster_disease_modules:
    input:
        network="results/networks/disease_{disease}.tsv"
    output:
        modules="results/modules/modules_diseasePPI_{disease}.tsv",
        meta="results/modules/modules_diseasePPI_{disease}.tsv.meta.json"
    params:
        mode="disease_combined"
    conda:
        "../../envs/clustering.yaml"
    script:
        "../scripts/cluster_graph.py"
