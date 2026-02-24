
import sys

try:
    from pathlib import Path

    sys.path.insert(0, str(Path(snakemake.scriptdir)))
except Exception:
    pass


import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from utils.io import read_tsv, write_tsv
from utils.provenance import write_metadata

cfg = snakemake.config
rule_name = snakemake.rule

in_leiden = snakemake.input["string_leiden"]
in_louvain = snakemake.input["string_louvain"]
in_disease = snakemake.input["disease"]
out_cmp = snakemake.output["comparison"]
out_meta = snakemake.output["meta"]


def _labels(df: pd.DataFrame, method: str) -> pd.Series:
    x = df[df["method"] == method].copy()
    x = x.drop_duplicates(subset=["gene_symbol"])
    return x.set_index("gene_symbol")["module_id"]


sl = read_tsv(in_leiden)
sv = read_tsv(in_louvain)
df = read_tsv(in_disease)

l1 = _labels(sl, "leiden")
l2 = _labels(sv, "louvain")
d1 = _labels(df, "leiden") if "leiden" in set(df.get("method", [])) else pd.Series(dtype=object)
d2 = _labels(df, "louvain") if "louvain" in set(df.get("method", [])) else pd.Series(dtype=object)

pairs = [
    ("string_leiden", l1, "disease_leiden", d1),
    ("string_leiden", l1, "disease_louvain", d2),
    ("string_louvain", l2, "disease_leiden", d1),
    ("string_louvain", l2, "disease_louvain", d2),
]

rows = []
for name_a, a, name_b, b in pairs:
    common = sorted(set(a.index) & set(b.index))
    if len(common) < 2:
        rows.append({"partition_a": name_a, "partition_b": name_b, "n_common": len(common), "ari": "", "nmi": ""})
        continue
    la = [a[g] for g in common]
    lb = [b[g] for g in common]
    rows.append(
        {
            "partition_a": name_a,
            "partition_b": name_b,
            "n_common": len(common),
            "ari": adjusted_rand_score(la, lb),
            "nmi": normalized_mutual_info_score(la, lb),
        }
    )

out = pd.DataFrame(rows)
write_tsv(out, out_cmp)

write_metadata(
    output_meta_path=out_meta,
    rule_name=rule_name,
    config=cfg,
    input_paths=[in_leiden, in_louvain, in_disease],
    parameters_used={"metrics": ["ARI", "NMI"]},
    organism_taxid=cfg["project"]["organism_taxid"],
    software_versions={"python": sys.version.split()[0], "pandas": pd.__version__},
)
