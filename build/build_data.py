#!/usr/bin/env python3
"""Build the source-traceable static release of the Human Aging Atlas."""

from __future__ import annotations

import json
import shutil
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from source_loaders import (
    ASSETS_DIR,
    CHUNK_SIZE,
    CHUNKS_DIR,
    DATA_DIR,
    DOWNLOADS,
    EPIGENETIC_DATA_URL,
    EPIGENETIC_PUBLICATION_URL,
    ORGANAGE_COMMIT,
    ORGANAGE_MODEL_ROOT_URL,
    ORGANAGE_PUBLICATION_URL,
    PUBLIC_ATLAS,
    SNAPSHOT_DIR,
    TRANSCRIPTOMIC_DATA_URL,
    TRANSCRIPTOMIC_PUBLICATION_URL,
    add_gene_evidence_cpgs,
    attach_reactome_metabolites,
    attach_reactome_proteins,
    clean,
    compact_dict,
    dump_json,
    evidence_layer_summary,
    gene_selection_metrics,
    load_chebi,
    load_hgnc,
    load_metaboage,
    load_public_genes,
    load_reactome,
    load_snapshot,
    load_uniprot,
    natural_key,
    prune_gene_record,
    rank_gene_candidates,
    sha256,
    source_definitions,
)


GENE_LIMIT = 1000
SCHEMA_VERSION = 10
PRIORITY_GENE_SYMBOLS = (
    "TP53",
    "APOE",
    "FOXO3",
    "FOXO1",
    "SIRT1",
    "SIRT3",
    "MTOR",
    "IGF1",
    "IGF1R",
    "KL",
    "TERT",
    "LMNA",
    "CDKN2A",
    "CDKN1A",
    "WRN",
    "ATM",
    "SOD2",
    "AKT1",
    "PPARGC1A",
    "NFE2L2",
    "RB1",
    "IL6",
    "TNF",
)
MODULES = (
    {
        "key": "genomics",
        "label": "Genomics",
        "primaryType": "gene",
        "sources": ["GenAge", "LongevityMap"],
        "description": "Curated aging genes and human longevity associations.",
    },
    {
        "key": "epigenomics",
        "label": "Epigenomics",
        "primaryType": "cpg",
        "sources": ["cAge", "bAge"],
        "description": "Gene-annotated CpG sites associated with age or mortality.",
    },
    {
        "key": "transcriptomics",
        "label": "Transcriptomics",
        "primaryType": "gene",
        "sources": ["tAge"],
        "description": "Human multi-tissue expression associations with age.",
    },
    {
        "key": "proteomics",
        "label": "Proteomics",
        "primaryType": "protein",
        "sources": ["OrganAge", "UniProtKB"],
        "description": "Human protein records with OrganAge evidence.",
    },
    {
        "key": "metabolomics",
        "label": "Metabolomics",
        "primaryType": "metabolite",
        "sources": ["MetaboAge DB", "ChEBI"],
        "description": "MetaboAge metabolites with source identity and directly mapped pathway context.",
    },
    {
        "key": "pathways",
        "label": "Pathways",
        "primaryType": "pathway",
        "sources": ["Reactome"],
        "description": "Direct human pathway mappings for retained atlas entities.",
    },
)


def choose_genes(
    genes: dict[str, dict[str, Any]],
    ranked: list[tuple[str, dict[str, Any]]],
) -> tuple[list[str], list[dict[str, Any]]]:
    ranked_symbols = [symbol for symbol, _ in ranked]
    priority = [
        symbol for symbol in PRIORITY_GENE_SYMBOLS if symbol in ranked_symbols
    ]
    selected = set(priority)
    for symbol in ranked_symbols:
        if len(selected) >= GENE_LIMIT:
            break
        selected.add(symbol)
    selected_symbols = [
        symbol for symbol in ranked_symbols if symbol in selected
    ]
    if len(selected_symbols) != GENE_LIMIT:
        raise RuntimeError(
            f"Expected {GENE_LIMIT} genes, found {len(selected_symbols)}"
        )
    audit = []
    priority_set = set(priority)
    for rank, symbol in enumerate(selected_symbols, 1):
        audit.append(
            {
                "releaseRank": rank,
                "symbol": symbol,
                "selectionBasis": (
                    "Priority aging and longevity reference gene"
                    if symbol in priority_set
                    else "Evidence-breadth ranking"
                ),
                **gene_selection_metrics(genes[symbol]),
            }
        )
    return selected_symbols, audit


def apply_gene_descriptions(
    genes: dict[str, dict[str, Any]],
    hgnc_by_symbol: dict[str, dict[str, str]],
) -> None:
    ncbi_snapshot = load_snapshot("ncbi_gene_descriptions.json")
    ncbi_records = ncbi_snapshot.get("records", {})
    for symbol, gene in genes.items():
        annotation = gene.get("annotation", {})
        entrez_id = clean(annotation.get("humanEntrezId"))
        ncbi = ncbi_records.get(entrez_id or "", {})
        summary = clean(ncbi.get("summary"))
        if summary:
            source_url = (
                f"https://www.ncbi.nlm.nih.gov/gene/{entrez_id}"
                if entrez_id
                else annotation.get("ncbiUrl")
            )
            gene["summary"] = summary
            gene["summarySource"] = {
                "label": "NCBI Gene summary",
                "url": source_url,
                "humanEntrezId": entrez_id,
            }
            continue

        source_name = (
            clean(ncbi.get("description"))
            or clean(hgnc_by_symbol.get(symbol, {}).get("name"))
            or clean(annotation.get("approvedName"))
        )
        if not source_name:
            raise RuntimeError(f"{symbol}: no source-backed gene description")
        gene["summary"] = source_name
        gene["summarySource"] = {
            "label": "HGNC approved gene description",
            "url": annotation.get("hgncUrl"),
        }


def add_evidence_urls(genes: dict[str, dict[str, Any]]) -> None:
    for symbol, gene in genes.items():
        evidence = gene.get("evidence", {})
        genage = evidence.get("genAgeHuman")
        if genage:
            genage["sourceUrl"] = (
                "https://genomics.senescence.info/genes/entry.php?"
                f"hgnc={symbol}"
            )
            genage.pop("geneName", None)
        for record in evidence.get("longevityMap", []):
            record["sourceUrl"] = (
                "https://genomics.senescence.info/longevity/entry.php?"
                f"id={record['reportId']}"
            )
        for key in ("epigeneticAge", "epigeneticMortality"):
            for record in evidence.get(key, []):
                record["sourceUrl"] = EPIGENETIC_DATA_URL
                record["publicationUrl"] = EPIGENETIC_PUBLICATION_URL
        for record in evidence.get("transcriptomic", []):
            record["sourceUrl"] = TRANSCRIPTOMIC_DATA_URL
            record["publicationUrl"] = TRANSCRIPTOMIC_PUBLICATION_URL
            source_symbol = record.pop("sourceMouseSymbol", None)
            if source_symbol:
                record["sourceSymbol"] = source_symbol
            for field in (
                "sourceMouseEntrezId",
                "orthologyClassId",
            ):
                record.pop(field, None)
        for record in evidence.get("organAge", []):
            if record.get("sourceCommit") != ORGANAGE_COMMIT:
                raise RuntimeError(
                    f"{symbol}: unexpected OrganAge source commit "
                    f"{record.get('sourceCommit')}"
                )
            organ = str(record.get("organ") or "").replace(" ", "%20")
            record["sourceUrl"] = f"{ORGANAGE_MODEL_ROOT_URL}/{organ}"
            record["publicationUrl"] = ORGANAGE_PUBLICATION_URL


def deduplicate_links(
    rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    return sorted(
        {row["id"]: row for row in rows}.values(),
        key=lambda row: (row.get("name", ""), natural_key(row["id"])),
    )


def field_coverage(records: dict[str, dict[str, Any]]) -> dict[str, Any]:
    fields = Counter(
        key
        for record in records.values()
        for key, value in record.items()
        if value not in (None, "", [], {})
    )
    total = len(records)
    return {
        field: {
            "present": count,
            "total": total,
            "fraction": round(count / total, 4) if total else 0,
        }
        for field, count in sorted(fields.items())
    }


def build() -> dict[str, Any]:
    if DATA_DIR.exists():
        shutil.rmtree(DATA_DIR)
    DATA_DIR.mkdir(parents=True)
    CHUNKS_DIR.mkdir(parents=True)
    ASSETS_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(
        PUBLIC_ATLAS / "assets" / "atlas-logo.png",
        ASSETS_DIR / "atlas-logo.png",
    )

    all_genes = load_public_genes()
    ranked = rank_gene_candidates(all_genes)
    selected_symbols, selection_audit = choose_genes(all_genes, ranked)
    selected = set(selected_symbols)
    hgnc_by_symbol, by_entrez = load_hgnc()
    genes = {
        symbol: prune_gene_record(all_genes[symbol])
        for symbol in selected_symbols
    }
    apply_gene_descriptions(genes, hgnc_by_symbol)
    add_evidence_urls(genes)

    pathways, gene_pathways = load_reactome(selected, by_entrez)
    proteins, gene_proteins = load_uniprot(selected, genes)
    proteins = {
        accession: protein
        for accession, protein in proteins.items()
        if genes[protein["geneSymbol"]]
        .get("evidence", {})
        .get("organAge")
        and protein.get("function")
    }
    gene_proteins = {
        symbol: [
            accession
            for accession in accessions
            if accession in proteins
        ]
        for symbol, accessions in gene_proteins.items()
    }
    attach_reactome_proteins(pathways, proteins)
    for protein in proteins.values():
        symbol = protein["geneSymbol"]
        protein["pathways"] = deduplicate_links(
            protein.get("pathways", [])
        )
        protein["organAgeEvidence"] = genes[symbol].get(
            "evidence", {}
        ).get("organAge", [])

    compounds, metabolite_name_lookup = load_chebi()
    attach_reactome_metabolites(
        pathways,
        compounds,
    )
    aging_metabolites = load_metaboage(
        compounds,
        metabolite_name_lookup,
    )
    metabolites = {
        entity_id: compounds[entity_id]
        for entity_id in sorted(aging_metabolites, key=natural_key)
        if entity_id in compounds
        and compounds[entity_id].get("definition")
    }
    for metabolite in metabolites.values():
        metabolite["pathways"] = deduplicate_links(
            [
                row
                for row in metabolite.get("pathways", [])
                if row["id"] in pathways
            ]
        )

    cpgs, gene_cpgs = add_gene_evidence_cpgs(genes)

    for pathway in pathways.values():
        pathway["genes"] = [
            {
                **member,
                "name": genes[member["symbol"]]
                .get("annotation", {})
                .get("approvedName"),
                "evidenceLayers": evidence_layer_summary(
                    genes[member["symbol"]]
                ),
            }
            for member in sorted(
                {
                    row["symbol"]: row
                    for row in pathway.get("genes", [])
                    if row["symbol"] in genes
                }.values(),
                key=lambda row: row["symbol"],
            )
        ]
        pathway["proteins"] = sorted(
            {
                row["id"]: row
                for row in pathway.get("proteins", [])
                if row["id"] in proteins
            }.values(),
            key=lambda row: natural_key(row["id"]),
        )
        pathway["metabolites"] = sorted(
            {
                row["id"]: row
                for row in pathway.get("metabolites", [])
                if row["id"] in metabolites
            }.values(),
            key=lambda row: (row.get("name", ""), natural_key(row["id"])),
        )

    gene_details: dict[str, dict[str, Any]] = {}
    for symbol in selected_symbols:
        gene = genes[symbol]
        protein_ids = sorted(
            {
                accession
                for accession in gene_proteins.get(symbol, [])
                if accession in proteins
            },
            key=natural_key,
        )
        pathway_links = deduplicate_links(
            [
                row
                for row in gene_pathways.get(symbol, [])
                if row["id"] in pathways
            ]
        )
        summary_source = gene["summarySource"]
        gene_details[symbol] = {
            **gene,
            "id": symbol,
            "entityType": "gene",
            "name": gene["annotation"]["approvedName"],
            "description": gene["summary"],
            "descriptionSource": summary_source["label"],
            "evidenceLayers": evidence_layer_summary(gene),
            "connections": {
                "proteins": protein_ids,
                "cpgs": sorted(
                    set(gene_cpgs.get(symbol, [])),
                    key=natural_key,
                ),
                "pathways": pathway_links,
            },
            "provenance": [
                {
                    "source": summary_source["label"],
                    "relationship": "Gene description",
                    "sourceUrl": summary_source["url"],
                },
                {
                    "source": "HGNC",
                    "relationship": "Approved human gene identity",
                    "sourceUrl": gene["annotation"]["hgncUrl"],
                },
            ],
        }

    entities: dict[str, dict[str, dict[str, Any]]] = {
        "gene": gene_details,
        "cpg": cpgs,
        "protein": proteins,
        "metabolite": metabolites,
        "pathway": pathways,
    }

    manifest_entities: dict[str, dict[str, Any]] = {}
    search_index: list[dict[str, Any]] = []
    for entity_type, records in entities.items():
        ordered_ids = sorted(records, key=natural_key)
        chunks: list[str] = []
        id_to_chunk: dict[str, str] = {}
        for index in range(0, len(ordered_ids), CHUNK_SIZE):
            chunk_ids = ordered_ids[index : index + CHUNK_SIZE]
            chunk_name = f"{entity_type}-{index // CHUNK_SIZE:03d}.json"
            dump_json(
                CHUNKS_DIR / chunk_name,
                {entity_id: records[entity_id] for entity_id in chunk_ids},
            )
            chunks.append(chunk_name)
            for entity_id in chunk_ids:
                id_to_chunk[entity_id] = chunk_name
        manifest_entities[entity_type] = {
            "count": len(records),
            "chunks": chunks,
            "idToChunk": id_to_chunk,
        }

        for entity_id, record in records.items():
            aliases: list[str] = []
            layers: list[dict[str, Any]] = []
            modules: list[str] = []
            details: dict[str, Any] = {}
            if entity_type == "gene":
                annotation = record["annotation"]
                aliases = (
                    annotation.get("aliases", [])
                    + annotation.get("previousSymbols", [])
                )
                layers = record["evidenceLayers"]
                modules = [layer["key"] for layer in layers]
                details = {
                    "location": annotation.get("chromosomeLocation"),
                    "sourceUrl": record["summarySource"]["url"],
                    "evidenceSourceUrls": {
                        **(
                            {
                                "GenAge": record["evidence"][
                                    "genAgeHuman"
                                ]["sourceUrl"]
                            }
                            if record["evidence"].get("genAgeHuman")
                            else {}
                        ),
                        **(
                            {
                                "LongevityMap": record["evidence"][
                                    "longevityMap"
                                ][0]["sourceUrl"]
                            }
                            if record["evidence"].get("longevityMap")
                            else {}
                        ),
                        **(
                            {
                                "tAge": record["evidence"][
                                    "transcriptomic"
                                ][0]["sourceUrl"]
                            }
                            if record["evidence"].get("transcriptomic")
                            else {}
                        ),
                    },
                }
            elif entity_type == "cpg":
                modules = ["epigenomics"]
                layers = [
                    {
                        "key": "epigenomics",
                        "sources": sorted(
                            {
                                row["source"]
                                for row in record["agingEvidence"]
                            }
                        ),
                    }
                ]
                aliases = [
                    row["symbol"] for row in record["geneMappings"]
                ]
                details = {
                    "geneSymbols": aliases,
                    "endpoints": sorted(
                        {
                            row["endpoint"]
                            for row in record["agingEvidence"]
                        }
                    ),
                    "sourceUrl": record["url"],
                }
            elif entity_type == "protein":
                modules = ["proteomics"]
                aliases = [
                    record.get("entryName", ""),
                    record.get("geneSymbol", ""),
                ]
                if record.get("organAgeEvidence"):
                    layers = [
                        {"key": "proteomics", "sources": ["OrganAge"]}
                    ]
                details = {
                    "geneSymbol": record.get("geneSymbol"),
                    "hasOrganAgeEvidence": bool(
                        record.get("organAgeEvidence")
                    ),
                    "hasPathways": bool(record.get("pathways")),
                    "sourceUrl": record["url"],
                }
            elif entity_type == "metabolite":
                modules = ["metabolomics"]
                aliases = record.get("synonyms", [])[:30]
                if record.get("agingEvidence"):
                    layers = [
                        {
                            "key": "metabolomics",
                            "sources": ["MetaboAge DB"],
                        }
                    ]
                details = {
                    "hasAgingEvidence": bool(record.get("agingEvidence")),
                    "hasPathways": bool(record.get("pathways")),
                    "sourceUrl": record["url"],
                }
            elif entity_type == "pathway":
                modules = ["pathways"]
                aliases = [
                    member["symbol"]
                    for member in record.get("genes", [])
                ]
                details = {
                    "memberTypes": [
                        label
                        for label, rows in (
                            ("Genes", record.get("genes", [])),
                            ("Proteins", record.get("proteins", [])),
                            ("Metabolites", record.get("metabolites", [])),
                        )
                        if rows
                    ],
                    "sourceUrl": record["recordUrl"],
                }

            search_index.append(
                compact_dict(
                    {
                        "id": entity_id,
                        "type": entity_type,
                        "name": record.get("name") or entity_id,
                        "description": record.get("description"),
                        "aliases": [item for item in aliases if item],
                        "evidenceLayers": layers,
                        "modules": modules,
                        "details": details,
                        "chunk": id_to_chunk[entity_id],
                    }
                )
            )

    collections = {
        "modules": list(MODULES),
        "sources": source_definitions(),
    }
    dump_json(DATA_DIR / "collections.json", collections, pretty=True)
    dump_json(DATA_DIR / "search-index.json", search_index)
    dump_json(DATA_DIR / "selection-audit.json", selection_audit, pretty=True)
    dump_json(
        DATA_DIR / "home-summary.json",
        {
            "entityCounts": {
                entity_type: len(records)
                for entity_type, records in entities.items()
            }
        },
        pretty=True,
    )
    dump_json(
        DATA_DIR / "quality-report.json",
        {
            "generatedAt": datetime.now(UTC).isoformat(timespec="seconds"),
            "fieldCoverage": {
                entity_type: field_coverage(records)
                for entity_type, records in entities.items()
            },
            "sourceDescriptionCoverage": {
                entity_type: {
                    "withDescription": sum(
                        bool(record.get("description"))
                        for record in records.values()
                    ),
                    "total": len(records),
                }
                for entity_type, records in entities.items()
                if entity_type
                in {"gene", "cpg", "protein", "metabolite", "pathway"}
            },
        },
        pretty=True,
    )

    source_paths = [
        PUBLIC_ATLAS / "data" / "build-report.json",
        DOWNLOADS / "01_GenAge" / "human_genes" / "genage_human.csv",
        DOWNLOADS / "02_LongevityMap" / "longevity_genes" / "longevity.csv",
        DOWNLOADS / "03_cAge" / "13073_2023_1161_MOESM4_ESM.xlsx",
        DOWNLOADS / "05_tAge" / "41586_2026_10542_MOESM4_ESM.xlsx",
        PUBLIC_ATLAS / "build" / "derived" / "organage_features.csv",
        DOWNLOADS / "13_MetaboAge" / "variation-data.xlsx",
        DOWNLOADS / "14_ChEBI" / "compounds.tsv.gz",
        DOWNLOADS / "14_ChEBI" / "names.tsv.gz",
        DOWNLOADS / "14_ChEBI" / "chemical_data.tsv.gz",
        DOWNLOADS / "16_Reactome" / "ReactomePathways.txt",
        DOWNLOADS / "16_Reactome" / "NCBI2Reactome_PE_Pathway.txt",
        DOWNLOADS / "16_Reactome" / "UniProt2Reactome_PE_Pathway.txt",
        DOWNLOADS / "16_Reactome" / "ChEBI2Reactome_PE_Pathway.txt",
        DOWNLOADS / "26_HGNC" / "hgnc_complete_set.txt",
        DOWNLOADS / "27_UniProt" / "uniprot_human_9606.tsv.gz",
        SNAPSHOT_DIR / "ncbi_gene_descriptions.json",
        SNAPSHOT_DIR / "reactome_pathways.json",
    ]
    manifest = {
        "schemaVersion": SCHEMA_VERSION,
        "generatedAt": datetime.now(UTC).isoformat(timespec="seconds"),
        "releaseType": "Reduced source-traceable internal static reference",
        "geneSelection": {
            "publishedCount": len(gene_details),
            "priorityPanel": list(PRIORITY_GENE_SYMBOLS),
            "fillRule": (
                "Descending retained evidence-layer breadth, source breadth, "
                "capped source-record depth, then gene symbol."
            ),
            "auditFile": "selection-audit.json",
        },
        "modules": list(MODULES),
        "entities": manifest_entities,
        "sourceChecksums": [
            {"path": str(path), "sha256": sha256(path)}
            for path in source_paths
            if path.exists()
        ],
        "scientificRules": [
            "Only human aging evidence is retained.",
            "Descriptions are copied or deterministically cleaned from cited sources; no generated scientific prose is published.",
            "Gene, protein, and metabolite pathway links require direct human Reactome physical-entity mappings.",
            "No gene-metabolite relationship is inferred from pathway co-membership.",
            "CpG gene annotation does not establish regulation or causality.",
            "Pathway context is not counted as an independent aging association.",
            "Empty sections are omitted and do not imply biological absence.",
            "No universal importance, causal, or clinical-actionability score is assigned.",
        ],
    }
    dump_json(DATA_DIR / "manifest.json", manifest, pretty=True)
    return {
        entity_type: len(records)
        for entity_type, records in entities.items()
    } | {"searchRecords": len(search_index)}


if __name__ == "__main__":
    print(json.dumps(build(), indent=2))
