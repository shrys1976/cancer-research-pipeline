
import hashlib
import json
import sys

try:
    from pathlib import Path

    sys.path.insert(0, str(Path(snakemake.scriptdir)))
except Exception:
    pass

import pandas as pd
import requests

from utils.io import ensure_parent, read_tsv, write_tsv
from utils.provenance import write_metadata

cfg = snakemake.config
rule_name = snakemake.rule

inputs = [snakemake.input["string_leiden"], snakemake.input["string_louvain"], snakemake.input["disease"]]
out_enrichment = snakemake.output["enrichment"]
out_meta = snakemake.output["meta"]

e_cfg = cfg["enrichment"]
api_url = e_cfg["api_url"]
organism = cfg["project"]["organism_gp"]
sources = e_cfg["sources"]
threshold = float(e_cfg["user_threshold"])
cache_dir = Path(e_cfg["cache_dir"])
cache_dir.mkdir(parents=True, exist_ok=True)


def _cache_key(payload: dict) -> str:
    raw = json.dumps(payload, sort_keys=True).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _query_gprofiler(genes: list[str]) -> list[dict]:
    payload = {
        "organism": organism,
        "query": genes,
        "sources": sources,
        "user_threshold": threshold,
    }
    key = _cache_key(payload)
    cpath = cache_dir / f"{key}.json"
    if cpath.exists():
        return json.loads(cpath.read_text(encoding="utf-8"))

    resp = requests.post(api_url, json=payload, timeout=60)
    resp.raise_for_status()
    data = resp.json().get("result", [])
    ensure_parent(cpath)
    cpath.write_text(json.dumps(data), encoding="utf-8")
    return data


rows = []
for path in inputs:
    df = read_tsv(path)
    if df.empty or "module_id" not in df.columns:
        continue

    methods = sorted(set(df["method"])) if "method" in df.columns else ["unknown"]
    for method in methods:
        sub = df[df["method"] == method] if "method" in df.columns else df
        for module_id, grp in sub.groupby("module_id"):
            genes = sorted(set(grp["gene_symbol"].astype(str)))
            if len(genes) < 2:
                continue
            try:
                terms = _query_gprofiler(genes)
            except Exception:
                terms = []
            for t in terms:
                rows.append(
                    {
                        "input_module_file": str(path),
                        "method": method,
                        "module_id": module_id,
                        "term_id": t.get("native", ""),
                        "term_name": t.get("name", ""),
                        "source": t.get("source", ""),
                        "p_value": t.get("p_value", ""),
                    }
                )

out = pd.DataFrame(rows)
if out.empty:
    out = pd.DataFrame(columns=["input_module_file", "method", "module_id", "term_id", "term_name", "source", "p_value"])

write_tsv(out, out_enrichment)

write_metadata(
    output_meta_path=out_meta,
    rule_name=rule_name,
    config=cfg,
    input_paths=inputs,
    parameters_used={"backend": "gprofiler", "organism": organism, "sources": sources, "user_threshold": threshold},
    organism_taxid=cfg["project"]["organism_taxid"],
    software_versions={"python": sys.version.split()[0], "pandas": pd.__version__, "requests": requests.__version__},
)
