
import json
import sys

try:
    from pathlib import Path

    sys.path.insert(0, str(Path(snakemake.scriptdir)))
except Exception:
    pass

import time

import pandas as pd
import requests

from utils.io import ensure_parent, read_tsv, write_tsv
from utils.provenance import write_metadata

cfg = snakemake.config
rule_name = snakemake.rule

in_seeds = snakemake.input["seeds"]
out_network = snakemake.output["network"]
out_meta = snakemake.output["meta"]

string_cfg = cfg["string"]
api_base = string_cfg["api_base"].rstrip("/")
fmt = string_cfg["format"]
score = int(string_cfg["required_score"])
partner_limit = int(string_cfg["partner_limit"])
timeout = int(string_cfg["timeout_seconds"])
retries = int(string_cfg["retries"])
cache_dir = Path(string_cfg["cache_dir"])
cache_dir.mkdir(parents=True, exist_ok=True)
species = cfg["project"]["organism_taxid"]


def _get_string_version() -> str:
    try:
        url = f"{api_base}/json/version"
        resp = requests.get(url, timeout=timeout)
        resp.raise_for_status()
        payload = resp.json()
        if isinstance(payload, list) and payload:
            maybe = payload[0].get("string_version") or payload[0].get("version")
            if maybe:
                return str(maybe)
        if isinstance(payload, dict):
            return str(payload.get("string_version") or payload.get("version") or "unknown")
        return "unknown"
    except Exception:  # noqa: BLE001
        return "unknown"


def _fetch_partners(identifier: str) -> pd.DataFrame:
    safe_id = identifier.replace("/", "_")
    cache_file = cache_dir / f"partners_{safe_id}_{species}_{score}.tsv"
    if cache_file.exists():
        return pd.read_csv(cache_file, sep="\t", dtype=str).fillna("")

    params = {
        "identifiers": identifier,
        "species": species,
        "required_score": score,
        "limit": partner_limit,
    }
    url = f"{api_base}/{fmt}/interaction_partners"

    last_err = None
    for attempt in range(retries):
        try:
            resp = requests.get(url, params=params, timeout=timeout)
            resp.raise_for_status()
            text = resp.text.strip()
            if not text:
                df = pd.DataFrame()
            else:
                from io import StringIO

                df = pd.read_csv(StringIO(text), sep="\t", dtype=str).fillna("")
            ensure_parent(cache_file)
            df.to_csv(cache_file, sep="\t", index=False)
            return df
        except Exception as err:  # noqa: BLE001
            last_err = err
            time.sleep(min(2**attempt, 5))
    raise RuntimeError(f"STRING API failed for {identifier}: {last_err}")


seeds = read_tsv(in_seeds)
protein_ids = sorted(set(seeds["protein_id"].tolist()))

edges = []
for pid in protein_ids:
    df = _fetch_partners(pid)
    if df.empty:
        continue

    for row in df.itertuples(index=False):
        node_a = getattr(row, "preferredName_A", "") or getattr(row, "stringId_A", "")
        node_b = getattr(row, "preferredName_B", "") or getattr(row, "stringId_B", "")
        s = getattr(row, "score", "")
        try:
            w = float(s)
        except Exception:  # noqa: BLE001
            w = 0.0
        if not node_a or not node_b:
            continue
        edges.append({"nodeA": node_a, "nodeB": node_b, "weight": w})

if not edges:
    raise RuntimeError("No STRING edges fetched. Check seed mapping and API availability.")

edf = pd.DataFrame(edges)
edf["key"] = edf.apply(lambda r: tuple(sorted((r["nodeA"], r["nodeB"]))), axis=1)
edf = edf.sort_values("weight", ascending=False).drop_duplicates("key")
edf[["nodeA", "nodeB"]] = pd.DataFrame(edf["key"].tolist(), index=edf.index)
edf = edf.drop(columns=["key"])
edf["source"] = "STRING"
string_version = _get_string_version()
edf["version"] = string_version
edf = edf.sort_values(["nodeA", "nodeB"]).reset_index(drop=True)

write_tsv(edf, out_network)

write_metadata(
    output_meta_path=out_meta,
    rule_name=rule_name,
    config=cfg,
    input_paths=[in_seeds],
    parameters_used={
        "required_score": score,
        "partner_limit": partner_limit,
        "string_version": string_version,
    },
    organism_taxid=species,
    software_versions={"python": sys.version.split()[0], "pandas": pd.__version__, "requests": requests.__version__},
)
