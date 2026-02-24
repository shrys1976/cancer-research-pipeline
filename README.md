# Cancer Research Pipeline (Snakemake)

Reproducible, config-driven disease-gene discovery pipeline implemented through Phase 8:

1. Setup check
2. Seed ingestion + UniProt canonical mapping
3. STRING background network retrieval
4. RWR expansion
5. Statistical proximity testing (degree-preserving null, 1000 perms, BH FDR)
6. Module detection (Leiden + Louvain)
7. Disease PPI retrieval from PSICQUIC
8. Disease PPI module detection
9. Module comparison (ARI, NMI)
10. Enrichment via g:Profiler API

## Run

```bash
snakemake -j 8 --use-conda
```

## Key Inputs

- Seed template: `resources/raw/cancer/cancer_seeds.tsv`
- Config: `config/config.yaml`
- Seed manifest: `manifests/cancer_seed_sources.yaml`

## Reproducibility

Each primary output emits a sidecar metadata JSON (`*.meta.json`) containing:

- `rule_name`
- `git_commit_hash`
- `config_hash`
- `input_files` / `input_versions`
- `parameters_used`
- `organism_taxid`
- `software_versions`
- `created_utc`
