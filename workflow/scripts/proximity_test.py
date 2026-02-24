
import math
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
in_expanded = snakemake.input["expanded"]
out_significant = snakemake.output["significant"]
out_meta = snakemake.output["meta"]

pr_cfg = cfg["proximity"]
n_perm = int(pr_cfg["n_permutations"])
seed_value = int(pr_cfg["random_seed"])
alpha = float(pr_cfg["alpha"])

seeds_df = read_tsv(in_seeds)
net_df = read_tsv(in_network)

G = nx.Graph()
for row in net_df.itertuples(index=False):
    G.add_edge(row.nodeA, row.nodeB, weight=float(row.weight))

seeds = sorted(set(seeds_df["gene_symbol"]))
seeds = [s for s in seeds if s in G]
if not seeds:
    raise RuntimeError("No seed nodes found in STRING graph for proximity test.")

obs = nx.multi_source_dijkstra_path_length(G, sources=seeds, weight=None)
obs_nodes = sorted(obs.keys())
obs_dist = {n: float(obs[n]) for n in obs_nodes}

count_le = {n: 0 for n in obs_nodes}

edges = G.number_of_edges()
nswap = max(edges, 1)

for i in range(n_perm):
    H = G.copy()
    try:
        nx.double_edge_swap(H, nswap=nswap, max_tries=nswap * 10, seed=seed_value + i)
    except Exception:
        pass
    rnd = nx.multi_source_dijkstra_path_length(H, sources=seeds, weight=None)
    for n in obs_nodes:
        rv = rnd.get(n, math.inf)
        if rv <= obs_dist[n]:
            count_le[n] += 1

rows = []
for n in obs_nodes:
    p = (count_le[n] + 1.0) / (n_perm + 1.0)
    rows.append({"gene_symbol": n, "observed_distance": obs_dist[n], "p_value": p})

res = pd.DataFrame(rows)
res = res[~res["gene_symbol"].isin(seeds)].copy()

res = res.sort_values("p_value").reset_index(drop=True)
m = len(res)
if m == 0:
    sig = pd.DataFrame(columns=["gene_symbol", "observed_distance", "p_value", "fdr_bh"])
else:
    rank = np.arange(1, m + 1)
    q = res["p_value"].to_numpy() * m / rank
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    res["fdr_bh"] = q
    sig = res[res["fdr_bh"] <= alpha].copy().sort_values(["fdr_bh", "p_value"])

write_tsv(sig, out_significant)

write_metadata(
    output_meta_path=out_meta,
    rule_name=rule_name,
    config=cfg,
    input_paths=[in_seeds, in_network, in_expanded],
    parameters_used={
        "null_model": "degree_preserving_double_edge_swap",
        "n_permutations": n_perm,
        "random_seed": seed_value,
        "fdr_method": "BH",
        "alpha": alpha,
    },
    organism_taxid=cfg["project"]["organism_taxid"],
    software_versions={"python": sys.version.split()[0], "numpy": np.__version__, "networkx": nx.__version__, "pandas": pd.__version__},
)
