rule compare_modules:
    input:
        string_leiden="results/modules/modules_{disease}_leiden.tsv",
        string_louvain="results/modules/modules_{disease}_louvain.tsv",
        disease="results/modules/modules_diseasePPI_{disease}.tsv"
    output:
        comparison="results/comparison/module_comparison_{disease}.tsv",
        meta="results/comparison/module_comparison_{disease}.tsv.meta.json"
    conda:
        "../../envs/comparison.yaml"
    script:
        "../scripts/compare_modules.py"
