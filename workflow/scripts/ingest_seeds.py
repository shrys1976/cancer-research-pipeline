
import sys

try:
    from pathlib import Path

    sys.path.insert(0, str(Path(snakemake.scriptdir)))
except Exception:
    pass

import pandas as pd

from utils.ids import map_gene_to_uniprot_canonical
from utils.io import read_tsv, write_tsv
from utils.provenance import write_metadata

cfg = snakemake.config
rule_name = snakemake.rule

in_seeds = snakemake.input["seeds_raw"]
out_seeds = snakemake.output["seeds"]
out_dropped = snakemake.output["dropped"]
out_meta = snakemake.output["meta"]

disease = snakemake.wildcards.disease
organism_taxid = cfg["project"]["organism_taxid"]
seed_cfg = cfg["seed_processing"]

df = read_tsv(in_seeds)
required_cols = ["gene_symbol", "protein_id", "weight", "source"]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing required seed columns: {missing}")

df["gene_symbol"] = df["gene_symbol"].str.strip()
df["source"] = df["source"].str.strip().replace("", "template")
df["weight"] = pd.to_numeric(df["weight"], errors="coerce").fillna(float(seed_cfg["default_weight"]))

mapped_rows = []
dropped_rows = []

for row in df.itertuples(index=False):
    gene_symbol = str(row.gene_symbol).strip()
    if not gene_symbol:
        dropped_rows.append({"gene_symbol": "", "reason": "empty_gene_symbol"})
        continue

    accession, status = map_gene_to_uniprot_canonical(
        gene_symbol=gene_symbol,
        organism_taxid=organism_taxid,
        uniprot_api=seed_cfg["uniprot_api"],
        timeout_seconds=int(seed_cfg["timeout_seconds"]),
        retries=int(seed_cfg["retries"]),
    )

    if accession is None:
        dropped_rows.append({"gene_symbol": gene_symbol, "reason": status})
        continue

    mapped_rows.append(
        {
            "gene_symbol": gene_symbol,
            "protein_id": accession,
            "weight": float(row.weight),
            "source": row.source or "template",
        }
    )

mapped = pd.DataFrame(mapped_rows)
if mapped.empty:
    raise RuntimeError("All seed genes dropped after mapping; check seed file and mapping constraints.")

mapped = mapped.drop_duplicates(subset=["gene_symbol", "protein_id"]).sort_values(["gene_symbol", "protein_id"])
write_tsv(mapped, out_seeds)

dropped = pd.DataFrame(dropped_rows)
if dropped.empty:
    dropped = pd.DataFrame(columns=["gene_symbol", "reason"])
write_tsv(dropped, out_dropped)

write_metadata(
    output_meta_path=out_meta,
    rule_name=rule_name,
    config=cfg,
    input_paths=[in_seeds, snakemake.input["manifest"]],
    parameters_used={
        "disease": disease,
        "mapping_policy": "reviewed_human_swissprot_canonical_only_drop_ambiguous",
    },
    organism_taxid=organism_taxid,
    software_versions={"python": sys.version.split()[0], "pandas": pd.__version__},
)
