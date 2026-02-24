
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd


def ensure_parent(path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def read_tsv(path: str | Path) -> "pd.DataFrame":
    import pandas as pd

    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def write_tsv(df: "pd.DataFrame", path: str | Path) -> None:
    ensure_parent(path)
    df.to_csv(path, sep="\t", index=False)
