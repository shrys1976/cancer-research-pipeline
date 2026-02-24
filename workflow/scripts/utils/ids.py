
import time
from typing import Any

import requests


def _request_with_retries(url: str, params: dict[str, Any], timeout: int, retries: int) -> dict:
    last_err = None
    for attempt in range(retries):
        try:
            resp = requests.get(url, params=params, timeout=timeout)
            resp.raise_for_status()
            return resp.json()
        except Exception as err:  # noqa: BLE001
            last_err = err
            time.sleep(min(2**attempt, 5))
    raise RuntimeError(f"UniProt request failed after {retries} retries: {last_err}")


def map_gene_to_uniprot_canonical(
    gene_symbol: str,
    organism_taxid: str,
    uniprot_api: str,
    timeout_seconds: int,
    retries: int,
) -> tuple[str | None, str]:
    params = {
        "query": f"(gene_exact:{gene_symbol}) AND (organism_id:{organism_taxid}) AND (reviewed:true)",
        "fields": "accession,gene_primary",
        "size": 10,
        "format": "json",
    }
    payload = _request_with_retries(uniprot_api, params, timeout_seconds, retries)
    results = payload.get("results", [])

    if not results:
        return None, "unresolved"

    accessions = []
    for entry in results:
        accession = entry.get("primaryAccession")
        if accession:
            accessions.append(accession)

    accessions = sorted(set(accessions))
    if len(accessions) != 1:
        return None, f"ambiguous_{len(accessions)}"

    return accessions[0], "mapped"
