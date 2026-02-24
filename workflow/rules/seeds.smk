rule ingest_seeds:
    input:
        seeds_raw=config["seed_processing"]["input_tsv"],
        manifest=config["project"]["seed_manifest"]
    output:
        seeds="results/genes/seeds_{disease}.tsv",
        dropped="results/reports/seeds_{disease}_dropped.tsv",
        meta="results/genes/seeds_{disease}.tsv.meta.json"
    conda:
        "../../envs/seeds.yaml"
    script:
        "../scripts/ingest_seeds.py"
