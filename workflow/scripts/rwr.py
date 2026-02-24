
import sys

try:
    from pathlib import Path

    sys.path.insert(0, str(Path(snakemake.scriptdir)))
except Exception:
    pass


import networkx as nx
import numpy as np
import pandas as pd

from utils.io import read_tsv, write_tsv
from utils.provenance import write_metadata

cfg = snakemake.config
rule_name = snakemake.rule

in_seeds = snakemake.input["seeds"]
in_network = snakemake.input["network"]
out_expanded = snakemake.output["expanded"]
out_meta = snakemake.output["meta"]

p_cfg = cfg["propagation"]
restart = float(p_cfg["restart_probability"])
max_iter = int(p_cfg["max_iter"])
tol = float(p_cfg["tol"])

seeds_df = read_tsv(in_seeds)
net_df = read_tsv(in_network)

G = nx.Graph()
for row in net_df.itertuples(index=False):
    w = float(row.weight)
    G.add_edge(row.nodeA, row.nodeB, weight=w)

seed_nodes = sorted(set(seeds_df["gene_symbol"]))
present_seeds = [s for s in seed_nodes if s in G]
if not present_seeds:
    raise RuntimeError("No seed gene symbols found in STRING graph nodes.")

nodes = list(G.nodes())
idx = {n: i for i, n in enumerate(nodes)}

A = nx.to_numpy_array(G, nodelist=nodes, weight="weight", dtype=float)
row_sums = A.sum(axis=1)
row_sums[row_sums == 0] = 1.0
W = A / row_sums[:, None]

p0 = np.zeros(len(nodes), dtype=float)
for s in present_seeds:
    p0[idx[s]] = 1.0
p0 /= p0.sum()

p = p0.copy()
for _ in range(max_iter):
    p_next = (1.0 - restart) * (W.T @ p) + restart * p0
    if np.linalg.norm(p_next - p, ord=1) < tol:
        p = p_next
        break
    p = p_next

out = pd.DataFrame(
    {
        "gene_symbol": nodes,
        "rwr_score": p,
        "is_seed": [1 if n in present_seeds else 0 for n in nodes],
        "method": "rwr",
    }
).sort_values("rwr_score", ascending=False)

write_tsv(out, out_expanded)

write_metadata(
    output_meta_path=out_meta,
    rule_name=rule_name,
    config=cfg,
    input_paths=[in_seeds, in_network],
    parameters_used={"restart_probability": restart, "max_iter": max_iter, "tol": tol},
    organism_taxid=cfg["project"]["organism_taxid"],
    software_versions={"python": sys.version.split()[0], "numpy": np.__version__, "networkx": nx.__version__, "pandas": pd.__version__},
)
