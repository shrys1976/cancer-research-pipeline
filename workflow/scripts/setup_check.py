
import platform
import sys
from pathlib import Path

try:
    sys.path.insert(0, str(Path(snakemake.scriptdir)))
except Exception:
    pass

from utils.provenance import write_metadata

rule_name = snakemake.rule
output_report = snakemake.output["report"]
output_meta = snakemake.output["meta"]
cfg = snakemake.config

html = f"""<html><body><h1>Setup Check</h1>
<p>Python: {sys.version}</p>
<p>Platform: {platform.platform()}</p>
<p>Disease: {cfg['project']['disease']}</p>
</body></html>"""

Path(output_report).parent.mkdir(parents=True, exist_ok=True)
with open(output_report, "w", encoding="utf-8") as fh:
    fh.write(html)

write_metadata(
    output_meta_path=output_meta,
    rule_name=rule_name,
    config=cfg,
    input_paths=[],
    parameters_used={"check": "environment_setup"},
    organism_taxid=cfg["project"]["organism_taxid"],
    software_versions={"python": sys.version.split()[0]},
)
