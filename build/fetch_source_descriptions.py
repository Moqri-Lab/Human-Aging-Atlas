#!/usr/bin/env python3
"""Snapshot exact NCBI Gene and Reactome descriptions used by the atlas."""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import UTC, datetime
from typing import Any, Iterable

from source_loaders import (
    SNAPSHOT_DIR,
    dump_json,
    load_hgnc,
    load_public_genes,
    load_reactome,
    rank_gene_candidates,
)


USER_AGENT = "Human-Aging-Atlas/3.0 (source snapshot)"
BATCH_SIZE = 150


def batches(values: list[str], size: int = BATCH_SIZE) -> Iterable[list[str]]:
    for index in range(0, len(values), size):
        yield values[index : index + size]


def request_bytes(
    url: str,
    *,
    data: bytes | None = None,
    content_type: str | None = None,
) -> bytes:
    headers = {"User-Agent": USER_AGENT}
    if content_type:
        headers["Content-Type"] = content_type
    request = urllib.request.Request(
        url,
        data=data,
        headers=headers,
        method="POST" if data is not None else "GET",
    )
    with urllib.request.urlopen(request, timeout=90) as response:
        return response.read()


def fetch_ncbi_gene_descriptions(gene_ids: list[str]) -> dict[str, Any]:
    records: dict[str, dict[str, str]] = {}
    endpoint = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    )
    for group in batches(gene_ids, 100):
        payload = urllib.parse.urlencode(
            {
                "db": "gene",
                "id": ",".join(group),
                "retmode": "xml",
                "tool": "human_aging_atlas",
            }
        ).encode()
        root = ET.fromstring(
            request_bytes(
                endpoint,
                data=payload,
                content_type="application/x-www-form-urlencoded",
            )
        )
        for gene in root.findall(".//Entrezgene"):
            gene_id = gene.findtext(".//Gene-track_geneid")
            if not gene_id:
                continue
            summary = gene.findtext("Entrezgene_summary")
            description = gene.findtext(".//Gene-ref_desc")
            symbol = gene.findtext(".//Gene-ref_locus")
            records[gene_id] = {
                key: value.strip()
                for key, value in {
                    "summary": summary,
                    "description": description,
                    "symbol": symbol,
                }.items()
                if value and value.strip()
            }
        time.sleep(0.35)
    return {
        "source": "NCBI Gene",
        "endpoint": endpoint,
        "retrievedAt": datetime.now(UTC).isoformat(timespec="seconds"),
        "records": records,
    }


def fetch_reactome_pathways(pathway_ids: list[str]) -> dict[str, Any]:
    records: dict[str, dict[str, str]] = {}
    endpoint = "https://reactome.org/ContentService/data/query/ids"
    for group in batches(pathway_ids, 20):
        payload = ",".join(group).encode()
        raw = request_bytes(
            endpoint,
            data=payload,
            content_type="text/plain",
        )
        response = json.loads(raw)
        for pathway in response:
            pathway_id = pathway.get("stId")
            if not pathway_id:
                continue
            summations = [
                item.get("text", "").strip()
                for item in pathway.get("summation", [])
                if item.get("text", "").strip()
            ]
            records[pathway_id] = {
                key: value
                for key, value in {
                    "name": pathway.get("displayName"),
                    "summary": " ".join(summations),
                    "schemaClass": pathway.get("schemaClass"),
                    "species": pathway.get("speciesName"),
                }.items()
                if value
            }
        time.sleep(0.2)
    return {
        "source": "Reactome ContentService",
        "endpoint": endpoint,
        "retrievedAt": datetime.now(UTC).isoformat(timespec="seconds"),
        "records": records,
    }


def main() -> None:
    genes = load_public_genes()
    selected_symbols = [
        symbol for symbol, _ in rank_gene_candidates(genes)[:1000]
    ]
    selected = set(selected_symbols)
    _, by_entrez = load_hgnc()
    pathways, _ = load_reactome(selected, by_entrez)
    gene_ids = sorted(
        {
            str(genes[symbol].get("annotation", {}).get("humanEntrezId"))
            for symbol in selected_symbols
            if genes[symbol].get("annotation", {}).get("humanEntrezId")
        },
        key=int,
    )
    pathway_ids = sorted(pathways)

    SNAPSHOT_DIR.mkdir(parents=True, exist_ok=True)
    ncbi_path = SNAPSHOT_DIR / "ncbi_gene_descriptions.json"
    existing_ncbi = (
        json.loads(ncbi_path.read_text(encoding="utf-8"))
        if ncbi_path.exists()
        else {}
    )
    ncbi = (
        existing_ncbi
        if len(existing_ncbi.get("records", {})) == len(gene_ids)
        else fetch_ncbi_gene_descriptions(gene_ids)
    )
    reactome = fetch_reactome_pathways(pathway_ids)
    dump_json(
        SNAPSHOT_DIR / "ncbi_gene_descriptions.json",
        ncbi,
        pretty=True,
    )
    dump_json(
        SNAPSHOT_DIR / "reactome_pathways.json",
        reactome,
        pretty=True,
    )
    print(
        json.dumps(
            {
                "selectedGenes": len(selected_symbols),
                "requestedNcbiGenes": len(gene_ids),
                "receivedNcbiGenes": len(ncbi["records"]),
                "requestedReactomePathways": len(pathway_ids),
                "receivedReactomePathways": len(reactome["records"]),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
