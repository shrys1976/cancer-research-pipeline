
import sys

try:
    from pathlib import Path

    sys.path.insert(0, str(Path(snakemake.scriptdir)))
except Exception:
    pass


import igraph as ig
import leidenalg
import networkx as nx
import pandas as pd
from community import community_louvain

from utils.io import read_tsv, write_tsv
from utils.provenance import write_metadata

cfg = snakemake.config
rule_name = snakemake.rule

mode = snakemake.params["mode"]
in_network = snakemake.input["network"]
out_modules = snakemake.output["modules"]
out_meta = snakemake.output["meta"]

seed = int(cfg["clustering"]["random_seed"])
resolution = float(cfg["clustering"]["resolution"])
methods = cfg["clustering"]["methods"]


def _graph_from_df(df: pd.DataFrame) -> nx.Graph:
    G = nx.Graph()
    for row in df.itertuples(index=False):
        G.add_edge(row.nodeA, row.nodeB, weight=float(row.weight))
    return G


def _cluster_method(G: nx.Graph, method: str) -> pd.DataFrame:
    if G.number_of_nodes() == 0:
        return pd.DataFrame(columns=["gene_symbol", "module_id", "method"])

    if method == "leiden":
        g = ig.Graph.TupleList(G.edges(), directed=False)
        part = leidenalg.find_partition(
            g,
            leidenalg.RBConfigurationVertexPartition,
            resolution_parameter=resolution,
            seed=seed,
        )
        rows = []
        for mod_id, members in enumerate(part):
            for v in members:
                rows.append({"gene_symbol": g.vs[v]["name"], "module_id": mod_id, "method": "leiden"})
        return pd.DataFrame(rows)

    if method == "louvain":
        partition = community_louvain.best_partition(G, weight="weight", random_state=seed, resolution=resolution)
        rows = [{"gene_symbol": n, "module_id": mid, "method": "louvain"} for n, mid in partition.items()]
        return pd.DataFrame(rows)

    raise ValueError(f"Unsupported method: {method}")


net_df = read_tsv(in_network)
G = _graph_from_df(net_df)

if mode == "string_significant":
    significant = read_tsv(snakemake.input["significant"])
    keep_nodes = set(significant["gene_symbol"]) & set(G.nodes())
    H = G.subgraph(keep_nodes).copy()
    method = snakemake.wildcards.method
    out = _cluster_method(H, method)
elif mode == "disease_combined":
    dfs = []
    for method in methods:
        dfs.append(_cluster_method(G, method))
    out = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame(columns=["gene_symbol", "module_id", "method"])
else:
    raise ValueError(f"Unknown cluster mode: {mode}")

out = out.sort_values(["method", "module_id", "gene_symbol"]).reset_index(drop=True)
write_tsv(out, out_modules)

inputs = [in_network]
if mode == "string_significant":
    inputs.append(snakemake.input["significant"])

write_metadata(
    output_meta_path=out_meta,
    rule_name=rule_name,
    config=cfg,
    input_paths=inputs,
    parameters_used={"mode": mode, "methods": methods, "resolution": resolution, "random_seed": seed},
    organism_taxid=cfg["project"]["organism_taxid"],
    software_versions={
        "python": sys.version.split()[0],
        "networkx": nx.__version__,
        "pandas": pd.__version__,
        "igraph": ig.__version__,
        "leidenalg": leidenalg.__version__,
    },
)
