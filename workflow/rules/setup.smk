rule setup_check:
    output:
        report="results/reports/setup_check.html",
        meta="results/reports/setup_check.html.meta.json"
    conda:
        "../../envs/setup.yaml"
    script:
        "../scripts/setup_check.py"
