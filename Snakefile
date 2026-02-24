configfile: "config/config.yaml"

DISEASE = config["project"]["disease"]
CLUSTER_METHODS = config["clustering"]["methods"]

include: "workflow/rules/setup.smk"
include: "workflow/rules/seeds.smk"
include: "workflow/rules/string.smk"
include: "workflow/rules/propagation.smk"
include: "workflow/rules/proximity.smk"
include: "workflow/rules/clustering.smk"
include: "workflow/rules/psicquic.smk"
include: "workflow/rules/comparison.smk"
include: "workflow/rules/enrichment.smk"

rule all:
    input:
        "results/reports/setup_check.html",
        f"results/genes/seeds_{DISEASE}.tsv",
        f"results/networks/string_bg_{DISEASE}.tsv",
        f"results/genes/expanded_{DISEASE}_rwr.tsv",
        f"results/genes/significant_{DISEASE}.tsv",
        expand(f"results/modules/modules_{DISEASE}_{{method}}.tsv", method=CLUSTER_METHODS),
        f"results/networks/disease_{DISEASE}.tsv",
        f"results/modules/modules_diseasePPI_{DISEASE}.tsv",
        f"results/comparison/module_comparison_{DISEASE}.tsv",
        f"results/enrichment/enrichment_{DISEASE}.tsv"
