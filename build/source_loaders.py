#!/usr/bin/env python3
"""Load the source records used by the reduced Human Aging Atlas release."""

from __future__ import annotations

import copy
import csv
import gzip
import hashlib
import html
import json
import math
import re
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parent.parent
WORKSPACE = PROJECT_DIR.parent
PUBLIC_ATLAS = WORKSPACE / "aging-evidence-atlas"
DOWNLOADS = WORKSPACE / "outputs" / "atlas_dataset_downloads" / "files"
DATA_DIR = PROJECT_DIR / "data"
CHUNKS_DIR = DATA_DIR / "chunks"
ASSETS_DIR = PROJECT_DIR / "assets"
SNAPSHOT_DIR = PROJECT_DIR / "source_snapshots"
CHUNK_SIZE = 25
EPIGENETIC_DATA_URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1186%2Fs13073-023-01161-y/MediaObjects/"
    "13073_2023_1161_MOESM4_ESM.xlsx"
)
EPIGENETIC_PUBLICATION_URL = (
    "https://doi.org/10.1186/s13073-023-01161-y"
)
TRANSCRIPTOMIC_DATA_URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs41586-026-10542-3/MediaObjects/"
    "41586_2026_10542_MOESM4_ESM.xlsx"
)
TRANSCRIPTOMIC_PUBLICATION_URL = (
    "https://doi.org/10.1038/s41586-026-10542-3"
)
ORGANAGE_COMMIT = "59303fd0dccc191be1ff34bf0bbf5efd8b90387a"
ORGANAGE_REPOSITORY_URL = (
    f"https://github.com/hamiltonoh/organage/tree/{ORGANAGE_COMMIT}"
)
ORGANAGE_MODEL_ROOT_URL = (
    f"{ORGANAGE_REPOSITORY_URL}/src/organage/data/ml_models/KADRC/"
    "Zprot_stableassayps_perf95"
)
ORGANAGE_PUBLICATION_URL = (
    "https://doi.org/10.1038/s41586-023-06802-1"
)
METABOAGE_DATA_URL = (
    "https://www.metaboage.info/static/website/variation-data.xlsx"
)
METABOAGE_DOWNLOAD_URL = "https://www.metaboage.info/download/"
METABOAGE_PUBLICATION_URL = (
    "https://doi.org/10.1007/s10522-020-09892-w"
)
RETAINED_EVIDENCE_KEYS = (
    "genAgeHuman",
    "longevityMap",
    "epigeneticAge",
    "epigeneticMortality",
    "transcriptomic",
    "organAge",
)


def clean(value: Any) -> str | None:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    return text or None


def plain_text(value: Any) -> str:
    text = html.unescape(clean(value) or "")
    text = re.sub(r"<[^>]+>", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def compact_dict(data: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in data.items()
        if value not in (None, "", [], {}, ())
    }


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if hasattr(value, "item") and not isinstance(value, (str, bytes)):
        try:
            value = value.item()
        except (TypeError, ValueError):
            pass
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def natural_key(value: str) -> tuple[Any, ...]:
    return tuple(
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", value)
    )


def normalized_name(value: Any) -> str:
    text = clean(value) or ""
    text = text.lower().replace("α", "alpha").replace("β", "beta")
    return re.sub(r"[^a-z0-9]+", "", text)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dump_json(path: Path, value: Any, *, pretty: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(
            json_safe(value),
            handle,
            ensure_ascii=True,
            indent=2 if pretty else None,
            separators=None if pretty else (",", ":"),
            allow_nan=False,
        )
        handle.write("\n")


def load_snapshot(name: str) -> dict[str, Any]:
    path = SNAPSHOT_DIR / name
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def load_public_genes() -> dict[str, dict[str, Any]]:
    genes: dict[str, dict[str, Any]] = {}
    for path in sorted((PUBLIC_ATLAS / "data").glob("genes-*.json")):
        genes.update(json.loads(path.read_text(encoding="utf-8")))
    return genes


def retained_source_flags(gene: dict[str, Any]) -> dict[str, bool]:
    evidence = gene.get("evidence", {})
    human_transcriptomic = [
        record
        for record in evidence.get("transcriptomic", [])
        if record.get("organism") == "Human"
    ]
    return {
        "GenAge": bool(evidence.get("genAgeHuman")),
        "LongevityMap": bool(evidence.get("longevityMap")),
        "cAge": bool(evidence.get("epigeneticAge")),
        "bAge": bool(evidence.get("epigeneticMortality")),
        "tAge": bool(human_transcriptomic),
        "OrganAge": bool(evidence.get("organAge")),
    }


def gene_selection_metrics(gene: dict[str, Any]) -> dict[str, Any]:
    evidence = gene.get("evidence", {})
    sources = retained_source_flags(gene)
    layers = {
        "genomics": sources["GenAge"] or sources["LongevityMap"],
        "epigenomics": sources["cAge"] or sources["bAge"],
        "transcriptomics": sources["tAge"],
        "proteomics": sources["OrganAge"],
    }
    record_count = 0
    for key in RETAINED_EVIDENCE_KEYS:
        value = evidence.get(key)
        if key == "transcriptomic":
            value = [
                record
                for record in value or []
                if record.get("organism") == "Human"
            ]
        if isinstance(value, dict):
            record_count += 1
        elif isinstance(value, list):
            record_count += len(value)
    return {
        "modalityCount": sum(layers.values()),
        "sourceCount": sum(sources.values()),
        "recordCount": record_count,
        "sources": [key for key, present in sources.items() if present],
        "modalities": [key for key, present in layers.items() if present],
    }


def rank_gene_candidates(
    genes: dict[str, dict[str, Any]],
) -> list[tuple[str, dict[str, Any]]]:
    ranked = []
    for symbol, gene in genes.items():
        annotation = gene.get("annotation", {})
        metrics = gene_selection_metrics(gene)
        if annotation.get("approvedName") and metrics["sourceCount"]:
            ranked.append((symbol, metrics))
    return sorted(
        ranked,
        key=lambda item: (
            -item[1]["modalityCount"],
            -item[1]["sourceCount"],
            -min(item[1]["recordCount"], 50),
            item[0],
        ),
    )


def prune_gene_record(gene: dict[str, Any]) -> dict[str, Any]:
    record = copy.deepcopy(gene)
    evidence = record.get("evidence", {})
    retained: dict[str, Any] = {}
    for key in RETAINED_EVIDENCE_KEYS:
        value = evidence.get(key)
        if key == "transcriptomic":
            value = [
                item
                for item in value or []
                if item.get("organism") == "Human"
            ]
        if value:
            retained[key] = value
    record["evidence"] = retained
    record.pop("mouseOrtholog", None)
    metrics = gene_selection_metrics(record)
    record["coverage"] = {
        "retainedSources": metrics["sources"],
        "retainedModalities": metrics["modalities"],
    }
    record.pop("statistics", None)
    return record


def load_hgnc() -> tuple[dict[str, dict[str, str]], dict[str, str]]:
    path = DOWNLOADS / "26_HGNC" / "hgnc_complete_set.txt"
    frame = pd.read_csv(path, sep="\t", dtype=str, low_memory=False).fillna("")
    frame = frame[frame["status"].eq("Approved")]
    by_symbol: dict[str, dict[str, str]] = {}
    by_entrez: dict[str, str] = {}
    for row in frame.to_dict(orient="records"):
        symbol = row["symbol"].upper()
        by_symbol[symbol] = row
        if row.get("entrez_id"):
            by_entrez[row["entrez_id"]] = symbol
    return by_symbol, by_entrez


def evidence_layer_summary(gene: dict[str, Any]) -> list[dict[str, Any]]:
    evidence = gene.get("evidence", {})
    layers: list[dict[str, Any]] = []
    genomics = []
    if evidence.get("genAgeHuman"):
        genomics.append("GenAge")
    if evidence.get("longevityMap"):
        genomics.append("LongevityMap")
    if genomics:
        layers.append({"key": "genomics", "sources": genomics})
    epigenomics = []
    if evidence.get("epigeneticAge"):
        epigenomics.append("cAge")
    if evidence.get("epigeneticMortality"):
        epigenomics.append("bAge")
    if epigenomics:
        layers.append({"key": "epigenomics", "sources": epigenomics})
    if any(
        row.get("organism") == "Human"
        for row in evidence.get("transcriptomic", [])
    ):
        layers.append({"key": "transcriptomics", "sources": ["tAge"]})
    if evidence.get("organAge"):
        layers.append({"key": "proteomics", "sources": ["OrganAge"]})
    return layers


def load_reactome(
    selected: set[str],
    by_entrez: dict[str, str],
) -> tuple[
    dict[str, dict[str, Any]],
    dict[str, list[dict[str, Any]]],
]:
    root = DOWNLOADS / "16_Reactome"
    descriptions = load_snapshot("reactome_pathways.json").get("records", {})
    names: dict[str, str] = {}
    with (root / "ReactomePathways.txt").open(encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 3 and row[2] == "Homo sapiens":
                names[row[0]] = row[1].strip()

    gene_members: dict[str, list[dict[str, Any]]] = defaultdict(list)
    gene_pathways: dict[str, set[str]] = defaultdict(set)
    mapping_path = root / "NCBI2Reactome_PE_Pathway.txt"
    with mapping_path.open(encoding="utf-8") as handle:
        for source_row, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(row) < 8 or row[7] != "Homo sapiens":
                continue
            symbol = by_entrez.get(row[0])
            pathway_id = row[3]
            if symbol not in selected or pathway_id not in names:
                continue
            entity_id = row[1]
            exact_url = (
                f"https://reactome.org/ContentService/data/query/{entity_id}"
            )
            pathway_url = (
                f"https://reactome.org/content/detail/{pathway_id}"
            )
            browser_url = f"{row[4]}&SEL={entity_id}"
            gene_pathways[symbol].add(pathway_id)
            gene_members[pathway_id].append(
                compact_dict(
                    {
                        "symbol": symbol,
                        "sourceIdentifier": row[0],
                        "physicalEntityId": entity_id,
                        "physicalEntity": row[2].strip(),
                        "evidenceCode": row[6],
                        "mappingScope": "Direct Reactome physical-entity mapping",
                        "sourceFile": mapping_path.name,
                        "sourceRow": source_row,
                        "sourceUrl": exact_url,
                        "pathwayUrl": pathway_url,
                        "browserUrl": browser_url,
                    }
                )
            )

    pathways: dict[str, dict[str, Any]] = {}
    for pathway_id in sorted(gene_members, key=natural_key):
        snapshot = descriptions.get(pathway_id, {})
        summary = plain_text(snapshot.get("summary"))
        name = plain_text(snapshot.get("name")) or names[pathway_id]
        description = summary or name
        pathways[pathway_id] = {
            "id": pathway_id,
            "entityType": "pathway",
            "name": name,
            "description": description,
            "descriptionSource": (
                "Reactome summation" if summary else "Reactome pathway name"
            ),
            "source": "Reactome",
            "species": "Homo sapiens",
            "url": (
                f"https://reactome.org/content/detail/{pathway_id}"
            ),
            "recordUrl": (
                f"https://reactome.org/ContentService/data/query/{pathway_id}"
            ),
            "genes": sorted(
                {
                    row["symbol"]: row
                    for row in gene_members[pathway_id]
                }.values(),
                key=lambda row: row["symbol"],
            ),
            "proteins": [],
            "metabolites": [],
            "provenance": [
                {
                    "source": "Reactome",
                    "relationship": "Direct human physical-entity membership",
                    "sourceUrl": (
                        "https://reactome.org/ContentService/data/query/"
                        f"{pathway_id}"
                    ),
                }
            ],
        }

    gene_links = {
        symbol: [
            {
                "id": pathway_id,
                "name": pathways[pathway_id]["name"],
                "source": "Reactome",
                **next(
                    member
                    for member in gene_members[pathway_id]
                    if member["symbol"] == symbol
                ),
            }
            for pathway_id in sorted(ids, key=lambda pid: names[pid])
            if pathway_id in pathways
        ]
        for symbol, ids in gene_pathways.items()
    }
    return pathways, gene_links


def parse_uniprot_function(value: Any) -> tuple[str | None, list[str]]:
    raw = clean(value)
    if not raw:
        return None, []
    pubmed_ids = sorted(
        set(re.findall(r"PubMed:(\d+)", raw)),
        key=natural_key,
    )
    text = re.sub(r"^FUNCTION:\s*", "", raw)
    text = re.sub(r"\s*\{[^{}]*\}\.?", "", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text or None, pubmed_ids


def recommended_protein_name(value: Any) -> str:
    text = clean(value) or ""
    match = re.match(r"^(.+?)(?=\s+\(|$)", text)
    return (match.group(1) if match else text).strip()


def load_uniprot(
    selected: set[str],
    genes: dict[str, dict[str, Any]],
) -> tuple[dict[str, dict[str, Any]], dict[str, list[str]]]:
    path = DOWNLOADS / "27_UniProt" / "uniprot_human_9606.tsv.gz"
    accession_to_symbol = {
        accession: symbol
        for symbol in selected
        for accession in genes[symbol]
        .get("annotation", {})
        .get("uniprotIds", [])
    }
    proteins: dict[str, dict[str, Any]] = {}
    gene_proteins: dict[str, list[str]] = defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = row.get("Entry", "")
            symbol = accession_to_symbol.get(accession)
            if not symbol:
                continue
            function, pubmed_ids = parse_uniprot_function(
                row.get("Function [CC]")
            )
            name = recommended_protein_name(row.get("Protein names"))
            description = function or name
            proteins[accession] = compact_dict(
                {
                    "id": accession,
                    "entityType": "protein",
                    "name": name,
                    "entryName": row.get("Entry Name"),
                    "geneSymbol": symbol,
                    "function": function,
                    "functionPubmedIds": pubmed_ids,
                    "description": description,
                    "descriptionSource": (
                        "UniProtKB Function annotation"
                        if function
                        else "UniProtKB recommended protein name"
                    ),
                    "url": (
                        f"https://www.uniprot.org/uniprotkb/{accession}/entry"
                    ),
                    "provenance": [
                        {
                            "source": "UniProtKB",
                            "relationship": "Protein identity and function",
                            "sourceUrl": (
                                "https://www.uniprot.org/uniprotkb/"
                                f"{accession}/entry"
                            ),
                        }
                    ],
                }
            )
            proteins[accession]["pathways"] = []
            proteins[accession]["organAgeEvidence"] = []
            gene_proteins[symbol].append(accession)
    return proteins, gene_proteins


def attach_reactome_proteins(
    pathways: dict[str, dict[str, Any]],
    proteins: dict[str, dict[str, Any]],
) -> None:
    path = DOWNLOADS / "16_Reactome" / "UniProt2Reactome_PE_Pathway.txt"
    with path.open(encoding="utf-8") as handle:
        for source_row, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(row) < 8 or row[7] != "Homo sapiens":
                continue
            accession, entity_id, physical_name, pathway_id = row[:4]
            if accession not in proteins or pathway_id not in pathways:
                continue
            source_url = (
                f"https://reactome.org/ContentService/data/query/{entity_id}"
            )
            pathway_url = (
                f"https://reactome.org/content/detail/{pathway_id}"
            )
            browser_url = f"{row[4]}&SEL={entity_id}"
            pathways[pathway_id]["proteins"].append(
                compact_dict(
                    {
                        "id": accession,
                        "name": physical_name,
                        "physicalEntityId": entity_id,
                        "evidenceCode": row[6],
                        "sourceFile": path.name,
                        "sourceRow": source_row,
                        "sourceUrl": source_url,
                        "pathwayUrl": pathway_url,
                        "browserUrl": browser_url,
                    }
                )
            )
            proteins[accession]["pathways"].append(
                {
                    "id": pathway_id,
                    "name": pathways[pathway_id]["name"],
                    "source": "Reactome",
                    "sourceUrl": source_url,
                    "pathwayUrl": pathway_url,
                    "browserUrl": browser_url,
                }
            )


def load_chebi() -> tuple[dict[str, dict[str, Any]], dict[str, str]]:
    root = DOWNLOADS / "14_ChEBI"
    compounds: dict[str, dict[str, Any]] = {}
    name_lookup: dict[str, str] = {}
    with gzip.open(root / "compounds.tsv.gz", "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            chebi_id = row.get("chebi_accession") or f"CHEBI:{row['id']}"
            name = plain_text(row.get("name"))
            definition = plain_text(row.get("definition"))
            compounds[chebi_id] = compact_dict(
                {
                    "id": chebi_id,
                    "entityType": "metabolite",
                    "name": name,
                    "definition": definition,
                    "description": definition or name,
                    "descriptionSource": (
                        "ChEBI definition"
                        if definition
                        else "ChEBI primary name"
                    ),
                    "url": (
                        "https://www.ebi.ac.uk/chebi/searchId.do?"
                        f"chebiId={chebi_id}"
                    ),
                    "provenance": [
                        {
                            "source": "ChEBI",
                            "relationship": "Metabolite identity",
                            "sourceUrl": (
                                "https://www.ebi.ac.uk/chebi/searchId.do?"
                                f"chebiId={chebi_id}"
                            ),
                        }
                    ],
                }
            )
            compounds[chebi_id]["pathways"] = []
            compounds[chebi_id]["agingEvidence"] = []
            compounds[chebi_id]["synonyms"] = []
            name_lookup[normalized_name(name)] = chebi_id
    with gzip.open(root / "names.tsv.gz", "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            chebi_id = f"CHEBI:{row.get('compound_id', '')}"
            name = plain_text(row.get("name"))
            if chebi_id in compounds and name:
                compounds[chebi_id]["synonyms"].append(name)
                name_lookup.setdefault(normalized_name(name), chebi_id)
    with gzip.open(
        root / "chemical_data.tsv.gz", "rt", encoding="utf-8"
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            chebi_id = f"CHEBI:{row.get('compound_id', '')}"
            if chebi_id not in compounds:
                continue
            compounds[chebi_id].update(
                compact_dict(
                    {
                        "formula": clean(row.get("formula")),
                        "charge": clean(row.get("charge")),
                        "averageMass": clean(row.get("mass")),
                        "monoisotopicMass": clean(
                            row.get("monoisotopic_mass")
                        ),
                    }
                )
            )
    return compounds, name_lookup


def attach_reactome_metabolites(
    pathways: dict[str, dict[str, Any]],
    compounds: dict[str, dict[str, Any]],
) -> set[str]:
    connected: set[str] = set()
    path = DOWNLOADS / "16_Reactome" / "ChEBI2Reactome_PE_Pathway.txt"
    with path.open(encoding="utf-8") as handle:
        for source_row, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(row) < 8 or row[7] != "Homo sapiens":
                continue
            raw_id, entity_id, physical_name, pathway_id = row[:4]
            if pathway_id not in pathways:
                continue
            chebi_id = (
                raw_id if raw_id.startswith("CHEBI:") else f"CHEBI:{raw_id}"
            )
            if chebi_id not in compounds:
                continue
            connected.add(chebi_id)
            name = re.sub(r"\s*\[[^\]]+\]\s*$", "", physical_name).strip()
            source_url = (
                f"https://reactome.org/ContentService/data/query/{entity_id}"
            )
            pathway_url = (
                f"https://reactome.org/content/detail/{pathway_id}"
            )
            browser_url = f"{row[4]}&SEL={entity_id}"
            pathways[pathway_id]["metabolites"].append(
                compact_dict(
                    {
                        "id": chebi_id,
                        "name": name,
                        "physicalEntityId": entity_id,
                        "evidenceCode": row[6],
                        "sourceFile": path.name,
                        "sourceRow": source_row,
                        "sourceUrl": source_url,
                        "pathwayUrl": pathway_url,
                        "browserUrl": browser_url,
                    }
                )
            )
            compounds[chebi_id]["pathways"].append(
                {
                    "id": pathway_id,
                    "name": pathways[pathway_id]["name"],
                    "source": "Reactome",
                    "sourceUrl": source_url,
                    "pathwayUrl": pathway_url,
                    "browserUrl": browser_url,
                }
            )
    return connected


def load_metaboage(
    compounds: dict[str, dict[str, Any]],
    name_lookup: dict[str, str],
) -> set[str]:
    root = DOWNLOADS / "13_MetaboAge"
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        chemistry = pd.read_excel(
            root / "chemical-modeling.xlsx",
            sheet_name="Sheet2",
        )
    variation = pd.read_excel(root / "variation-data.xlsx")
    chemistry_by_name = {
        normalized_name(row.get("Metabolite name")): row.to_dict()
        for _, row in chemistry.iterrows()
        if clean(row.get("Metabolite name"))
    }
    connected: set[str] = set()
    for source_row, row in variation.iterrows():
        name = clean(row.get("Metabolite"))
        if not name:
            continue
        key = normalized_name(name)
        chebi_id = name_lookup.get(key)
        if not chebi_id:
            chemistry_row = chemistry_by_name.get(key, {})
            raw_chebi = clean(
                chemistry_row.get("ChEBI ID")
                or chemistry_row.get("ChEBI")
                or chemistry_row.get("CHEBI ID")
            )
            match = re.search(r"(\d+)", raw_chebi or "")
            chebi_id = f"CHEBI:{match.group(1)}" if match else None
        if not chebi_id or chebi_id not in compounds:
            continue
        raw_pmid = clean(row.get("Study  PMID reference"))
        pmid_match = re.search(r"\d+", raw_pmid or "")
        publication_url = (
            f"https://pubmed.ncbi.nlm.nih.gov/{pmid_match.group(0)}/"
            if pmid_match
            else None
        )
        connected.add(chebi_id)
        compounds[chebi_id]["agingEvidence"].append(
            compact_dict(
                {
                    "source": "MetaboAge DB",
                    "metaboliteName": name,
                    "method": row.get("Method"),
                    "value": row.get("Mean value/beta value"),
                    "uncertainty": row.get(
                        "Standard deviation (+-)/standard error"
                    ),
                    "unit": row.get("Unit of measurment"),
                    "ageGroup": row.get("(NEW) Age group"),
                    "gender": row.get("Gender"),
                    "sample": row.get("Tipe of sample used in study"),
                    "pmid": pmid_match.group(0) if pmid_match else None,
                    "sourceFile": "variation-data.xlsx",
                    "sourceSheet": "Sheet1",
                    "sourceRow": int(source_row) + 2,
                    "sourceUrl": METABOAGE_DATA_URL,
                    "publicationUrl": publication_url,
                    "databaseUrl": METABOAGE_DOWNLOAD_URL,
                }
            )
        )
    return connected


def add_gene_evidence_cpgs(
    genes: dict[str, dict[str, Any]],
) -> tuple[dict[str, dict[str, Any]], dict[str, list[str]]]:
    cpgs: dict[str, dict[str, Any]] = {}
    gene_cpgs: dict[str, list[str]] = defaultdict(list)
    source_url = EPIGENETIC_DATA_URL
    source_workbook = (
        DOWNLOADS / "03_cAge" / "13073_2023_1161_MOESM4_ESM.xlsx"
    )
    source_titles = {
        sheet: plain_text(
            pd.read_excel(
                source_workbook,
                sheet_name=sheet,
                header=None,
                nrows=1,
            ).iat[0, 0]
        )
        for sheet in ("S1", "S3")
    }
    for symbol, gene in genes.items():
        evidence = gene.get("evidence", {})
        for key, source, endpoint in (
            ("epigeneticAge", "cAge", "Chronological age"),
            ("epigeneticMortality", "bAge", "All-cause mortality"),
        ):
            for record in evidence.get(key, []):
                cpg_id = clean(record.get("cpg"))
                if not cpg_id:
                    continue
                cpg = cpgs.setdefault(
                    cpg_id,
                    {
                        "id": cpg_id,
                        "entityType": "cpg",
                        "name": cpg_id,
                        "chromosome": clean(record.get("cpgChromosome")),
                        "position": record.get("cpgPosition"),
                        "genomeBuild": "GRCh38",
                        "geneMappings": [],
                        "agingEvidence": [],
                        "sourceDescriptions": [],
                        "url": source_url,
                        "provenance": [],
                    },
                )
                source_sheet = clean(record.get("sourceSheet"))
                source_description = source_titles.get(source_sheet or "")
                if source_description and not any(
                    item["source"] == source
                    for item in cpg["sourceDescriptions"]
                ):
                    cpg["sourceDescriptions"].append(
                        {
                            "source": source,
                            "sourceSheet": source_sheet,
                            "text": source_description,
                            "sourceUrl": source_url,
                            "publicationUrl": EPIGENETIC_PUBLICATION_URL,
                        }
                    )
                if not any(
                    mapping["symbol"] == symbol
                    for mapping in cpg["geneMappings"]
                ):
                    cpg["geneMappings"].append(
                        {
                            "symbol": symbol,
                            "mapping": "Source gene annotation",
                            "sourceUrl": source_url,
                        }
                    )
                cpg["agingEvidence"].append(
                    compact_dict(
                        {
                            "recordId": record.get("recordId"),
                            "source": source,
                            "endpoint": endpoint,
                            "model": record.get("model"),
                            "effectType": (
                                "Beta coefficient"
                                if key == "epigeneticAge"
                                else "Hazard ratio"
                            ),
                            "effect": (
                                record.get("beta")
                                if key == "epigeneticAge"
                                else record.get("hazardRatio")
                            ),
                            "standardError": record.get("standardError"),
                            "hazardRatioCiLow": record.get(
                                "hazardRatioCiLow"
                            ),
                            "hazardRatioCiHigh": record.get(
                                "hazardRatioCiHigh"
                            ),
                            "pValue": record.get("pValue"),
                            "sensitivityAnalysis": record.get(
                                "sensitivityAnalysis"
                            ),
                            "sourceSheet": record.get("sourceSheet"),
                            "sourceRow": record.get("sourceRow"),
                            "sourceUrl": source_url,
                            "publicationUrl": EPIGENETIC_PUBLICATION_URL,
                        }
                    )
                )
                gene_cpgs[symbol].append(cpg_id)
                if not any(
                    item["source"] == source for item in cpg["provenance"]
                ):
                    cpg["provenance"].append(
                        {
                            "source": source,
                            "relationship": endpoint,
                            "sourceUrl": source_url,
                        }
                    )
    for cpg in cpgs.values():
        description = cpg["sourceDescriptions"][0]
        cpg["description"] = description["text"]
        cpg["descriptionSource"] = (
            f"{description['sourceSheet']} source table title"
        )
    return cpgs, gene_cpgs


def source_definitions() -> list[dict[str, Any]]:
    return [
        {
            "name": "GenAge",
            "group": "Aging evidence",
            "layer": "genomics",
            "role": "Expert-curated aging genes",
            "url": "https://genomics.senescence.info/genes/",
        },
        {
            "name": "LongevityMap",
            "group": "Aging evidence",
            "layer": "genomics",
            "role": "Human longevity association reports",
            "url": "https://genomics.senescence.info/longevity/",
        },
        {
            "name": "cAge and bAge",
            "group": "Aging evidence",
            "layer": "epigenomics",
            "role": "Age- and mortality-associated CpG sites",
            "url": EPIGENETIC_DATA_URL,
            "publicationUrl": EPIGENETIC_PUBLICATION_URL,
        },
        {
            "name": "tAge",
            "group": "Aging evidence",
            "layer": "transcriptomics",
            "role": "Human multi-tissue transcriptomic associations",
            "url": TRANSCRIPTOMIC_DATA_URL,
            "publicationUrl": TRANSCRIPTOMIC_PUBLICATION_URL,
        },
        {
            "name": "OrganAge",
            "group": "Aging evidence",
            "layer": "proteomics",
            "role": "Organ-specific proteomic age models",
            "url": ORGANAGE_REPOSITORY_URL,
            "publicationUrl": ORGANAGE_PUBLICATION_URL,
        },
        {
            "name": "MetaboAge DB",
            "group": "Aging evidence",
            "layer": "metabolomics",
            "role": "Reported age-associated metabolite measurements",
            "url": METABOAGE_DATA_URL,
            "publicationUrl": METABOAGE_PUBLICATION_URL,
        },
        {
            "name": "Reactome",
            "group": "Connected context",
            "layer": "pathways",
            "role": "Direct human gene, protein, and metabolite pathway mappings",
            "url": "https://reactome.org/",
        },
        {
            "name": "NCBI Gene",
            "group": "Identity and description",
            "role": "Human gene summaries and identifiers",
            "url": "https://www.ncbi.nlm.nih.gov/gene/",
        },
        {
            "name": "HGNC",
            "group": "Identity and description",
            "role": "Approved human gene identifiers and names",
            "url": "https://www.genenames.org/",
        },
        {
            "name": "UniProtKB",
            "group": "Identity and description",
            "role": "Protein names, functions, and gene mappings",
            "url": "https://www.uniprot.org/",
        },
        {
            "name": "ChEBI",
            "group": "Identity and description",
            "role": "Metabolite names, definitions, and properties",
            "url": "https://www.ebi.ac.uk/chebi",
        },
    ]
