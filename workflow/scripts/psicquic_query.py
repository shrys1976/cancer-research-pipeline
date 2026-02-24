
import re
import sys

try:
    from pathlib import Path

    sys.path.insert(0, str(Path(snakemake.scriptdir)))
except Exception:
    pass

import pandas as pd
import requests

from utils.io import ensure_parent, write_tsv
from utils.provenance import write_metadata

cfg = snakemake.config
rule_name = snakemake.rule

out_network = snakemake.output["network"]
out_meta = snakemake.output["meta"]

p_cfg = cfg["psicquic"]
endpoint = p_cfg["endpoint"].rstrip("/")
queries = p_cfg["queries"]
max_lines = int(p_cfg["max_lines"])
cache_dir = Path(p_cfg["cache_dir"])
cache_dir.mkdir(parents=True, exist_ok=True)


def _extract_uniprot(token: str) -> str:
    for part in token.split("|"):
        m = re.match(r"uniprotkb:([A-Z0-9]+)", part)
        if m:
            return m.group(1)
    return ""


all_edges = []
version = "unknown"

for q in queries:
    query = f"species:9606 AND ({q})"
    cache_path = cache_dir / f"{q}.mitab"
    if cache_path.exists():
        text = cache_path.read_text(encoding="utf-8")
    else:
        url = f"{endpoint}/{query}?format=tab25"
        resp = requests.get(url, timeout=45)
        resp.raise_for_status()
        text = resp.text
        ensure_parent(cache_path)
        cache_path.write_text(text, encoding="utf-8")
        if "X-PSICQUIC-Spec-Version" in resp.headers:
            version = resp.headers.get("X-PSICQUIC-Spec-Version", "unknown")

    lines = [ln for ln in text.splitlines() if ln.strip()]
    for ln in lines[:max_lines]:
        cols = ln.split("\t")
        if len(cols) < 2:
            continue
        a = _extract_uniprot(cols[0])
        b = _extract_uniprot(cols[1])
        if not a or not b:
            continue
        all_edges.append({"nodeA": a, "nodeB": b, "weight": 1.0, "source": "PSICQUIC:IntAct", "version": version, "query": q})

if all_edges:
    df = pd.DataFrame(all_edges)
    df["key"] = df.apply(lambda r: tuple(sorted((r["nodeA"], r["nodeB"]))), axis=1)
    df = df.drop_duplicates("key").drop(columns=["key"])
else:
    df = pd.DataFrame(columns=["nodeA", "nodeB", "weight", "source", "version", "query"])

write_tsv(df.sort_values(["nodeA", "nodeB"]) if not df.empty else df, out_network)

write_metadata(
    output_meta_path=out_meta,
    rule_name=rule_name,
    config=cfg,
    input_paths=[],
    parameters_used={"queries": queries, "max_lines": max_lines, "psicquic_version": version},
    organism_taxid=cfg["project"]["organism_taxid"],
    software_versions={"python": sys.version.split()[0], "pandas": pd.__version__, "requests": requests.__version__},
)
