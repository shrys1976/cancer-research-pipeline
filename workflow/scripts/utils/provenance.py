
import datetime as dt
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

def ensure_parent(path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def _safe_git_commit() -> str:
    try:
        out = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
        return out
    except Exception:
        return "unknown"


def _sha256_jsonable(obj: Any) -> str:
    raw = json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _file_digest(path: str | Path) -> str:
    p = Path(path)
    if not p.exists():
        return "missing"
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_metadata(
    *,
    output_meta_path: str | Path,
    rule_name: str,
    config: dict,
    input_paths: list[str],
    parameters_used: dict,
    organism_taxid: str,
    software_versions: dict,
) -> None:
    meta = {
        "rule_name": rule_name,
        "git_commit_hash": _safe_git_commit(),
        "config_hash": _sha256_jsonable(config),
        "input_files": [{"path": p, "sha256": _file_digest(p)} for p in input_paths],
        "input_versions": {
            p: dt.datetime.utcfromtimestamp(Path(p).stat().st_mtime).isoformat() + "Z"
            for p in input_paths
            if Path(p).exists()
        },
        "parameters_used": parameters_used,
        "organism_taxid": str(organism_taxid),
        "software_versions": software_versions,
        "created_utc": dt.datetime.utcnow().isoformat() + "Z",
    }
    ensure_parent(output_meta_path)
    Path(output_meta_path).write_text(json.dumps(meta, indent=2), encoding="utf-8")
