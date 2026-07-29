#!/usr/bin/env python3
"""Validate the reduced source-traceable Human Aging Atlas release."""

from __future__ import annotations

import csv
import gzip
import html
import json
import math
import re
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parent.parent
WORKSPACE = PROJECT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
DOWNLOADS = WORKSPACE / "outputs" / "atlas_dataset_downloads" / "files"
PUBLIC_ATLAS = WORKSPACE / "aging-evidence-atlas"
ORGANAGE_COMMIT = "59303fd0dccc191be1ff34bf0bbf5efd8b90387a"
EXPECTED_TYPES = {"gene", "cpg", "protein", "metabolite", "pathway"}
EXPECTED_MODULES = {
    "genomics",
    "epigenomics",
    "transcriptomics",
    "proteomics",
    "metabolomics",
    "pathways",
}
FORBIDDEN_GENERATED_PHRASES = (
    "HGNC classifies",
    "Encoding-gene function summary",
    "Variant not reported",
)


def fail(message: str) -> None:
    raise AssertionError(message)


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def require_fields(record: dict[str, Any], fields: tuple[str, ...], context: str) -> None:
    missing = [
        field
        for field in fields
        if record.get(field) in (None, "", [], {})
    ]
    if missing:
        fail(f"{context}: missing {', '.join(missing)}")


def require_url(url: Any, context: str) -> None:
    parsed = urlparse(str(url or ""))
    if parsed.scheme != "https" or not parsed.netloc:
        fail(f"{context}: invalid source URL {url!r}")


def walk(value: Any):
    if isinstance(value, dict):
        yield value
        for child in value.values():
            yield from walk(child)
    elif isinstance(value, list):
        for child in value:
            yield from walk(child)


def source_scalar(value: Any) -> Any:
    if isinstance(value, dict) and "value" in value:
        value = value["value"]
    try:
        return None if pd.isna(value) else value
    except (TypeError, ValueError):
        return value


def normalized_text(value: Any) -> str:
    value = source_scalar(value)
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        value = int(value)
    return str(value).strip()


def plain_text(value: Any) -> str:
    text = html.unescape(normalized_text(value))
    text = re.sub(r"<[^>]+>", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def require_source_value(
    actual: Any,
    expected: Any,
    context: str,
) -> None:
    if isinstance(actual, dict):
        expected_text = normalized_text(expected)
        display = normalized_text(actual.get("display"))
        qualifier = normalized_text(actual.get("qualifier"))
        if display == expected_text:
            return
        if qualifier == "less_than" and expected_text.startswith("<"):
            require_source_value(
                actual.get("value"),
                expected_text.removeprefix("<"),
                context,
            )
            return
    actual = source_scalar(actual)
    expected = source_scalar(expected)
    if actual is None and expected is None:
        return
    try:
        actual_number = float(actual)
        expected_number = float(expected)
    except (TypeError, ValueError):
        if normalized_text(actual) != normalized_text(expected):
            fail(
                f"{context}: published {actual!r} differs from source "
                f"{expected!r}"
            )
        return
    if not math.isclose(
        actual_number,
        expected_number,
        rel_tol=1e-10,
        abs_tol=0.0,
    ):
        fail(
            f"{context}: published {actual!r} differs from source "
            f"{expected!r}"
        )


def source_row(
    frame: pd.DataFrame,
    row_number: Any,
    first_data_row: int,
    context: str,
) -> pd.Series:
    index = int(row_number) - first_data_row
    if index < 0 or index >= len(frame):
        fail(f"{context}: source row {row_number} is out of range")
    return frame.iloc[index]


def read_tsv_rows(path: Path) -> dict[int, list[str]]:
    with path.open(encoding="utf-8") as handle:
        return {
            row_number: row
            for row_number, row in enumerate(
                csv.reader(handle, delimiter="\t"),
                1,
            )
        }


def uniprot_function(value: Any) -> str:
    text = normalized_text(value)
    text = re.sub(r"^FUNCTION:\s*", "", text)
    text = re.sub(r"\s*\{[^{}]*\}\.?", "", text)
    return re.sub(r"\s+", " ", text).strip()


def validate_source_replay(
    entities: dict[str, dict[str, dict[str, Any]]],
) -> None:
    genage = pd.read_csv(
        DOWNLOADS / "01_GenAge" / "human_genes" / "genage_human.csv"
    )
    longevity = pd.read_csv(
        DOWNLOADS
        / "02_LongevityMap"
        / "longevity_genes"
        / "longevity.csv"
    )
    epigenetic_path = (
        DOWNLOADS / "03_cAge" / "13073_2023_1161_MOESM4_ESM.xlsx"
    )
    epigenetic = {
        sheet: pd.read_excel(
            epigenetic_path,
            sheet_name=sheet,
            header=2,
        )
        for sheet in ("S1", "S3")
    }
    epigenetic_titles = {
        sheet: plain_text(
            pd.read_excel(
                epigenetic_path,
                sheet_name=sheet,
                header=None,
                nrows=1,
            ).iat[0, 0]
        )
        for sheet in ("S1", "S3")
    }
    transcriptomic = pd.read_excel(
        DOWNLOADS / "05_tAge" / "41586_2026_10542_MOESM4_ESM.xlsx",
        sheet_name="(R) Human aging multi-tissue",
    )
    organage = pd.read_csv(
        PUBLIC_ATLAS / "build" / "derived" / "organage_features.csv"
    )
    metaboage = pd.read_excel(
        DOWNLOADS / "13_MetaboAge" / "variation-data.xlsx"
    )

    reactome_root = DOWNLOADS / "16_Reactome"
    reactome_rows = {
        name: read_tsv_rows(reactome_root / name)
        for name in (
            "NCBI2Reactome_PE_Pathway.txt",
            "UniProt2Reactome_PE_Pathway.txt",
            "ChEBI2Reactome_PE_Pathway.txt",
        )
    }
    reactome_snapshot = load_json(
        PROJECT_DIR / "source_snapshots" / "reactome_pathways.json"
    ).get("records", {})

    protein_ids = set(entities["protein"])
    uniprot_rows: dict[str, dict[str, str]] = {}
    with gzip.open(
        DOWNLOADS / "27_UniProt" / "uniprot_human_9606.tsv.gz",
        "rt",
        encoding="utf-8",
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("Entry") in protein_ids:
                uniprot_rows[row["Entry"]] = row

    metabolite_ids = set(entities["metabolite"])
    chebi_rows: dict[str, dict[str, str]] = {}
    with gzip.open(
        DOWNLOADS / "14_ChEBI" / "compounds.tsv.gz",
        "rt",
        encoding="utf-8",
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            chebi_id = row.get("chebi_accession") or f"CHEBI:{row['id']}"
            if chebi_id in metabolite_ids:
                chebi_rows[chebi_id] = row

    for symbol, gene in entities["gene"].items():
        evidence = gene.get("evidence", {})
        genage_record = evidence.get("genAgeHuman")
        if genage_record:
            row = source_row(
                genage,
                genage_record["sourceRow"],
                2,
                f"{symbol} GenAge",
            )
            require_source_value(
                genage_record.get("sourceSymbol"),
                row.get("symbol"),
                f"{symbol} GenAge symbol",
            )
            require_source_value(
                genage_record.get("genAgeId"),
                row.get("GenAge ID"),
                f"{symbol} GenAge ID",
            )

        for record in evidence.get("longevityMap", []):
            row = source_row(
                longevity,
                record["sourceRow"],
                2,
                f"{symbol} LongevityMap",
            )
            for field, source_field in (
                ("reportId", "id"),
                ("population", "Population"),
                ("variants", "Variant(s)"),
                ("pubmedId", "PubMed"),
            ):
                require_source_value(
                    record.get(field),
                    row.get(source_field),
                    f"{symbol} LongevityMap {field}",
                )

        for key in ("epigeneticAge", "epigeneticMortality"):
            for record in evidence.get(key, []):
                row = source_row(
                    epigenetic[record["sourceSheet"]],
                    record["sourceRow"],
                    4,
                    f"{symbol} {key}",
                )
                require_source_value(
                    record.get("cpg"),
                    row.get("CpG"),
                    f"{symbol} {key} CpG",
                )
                effect_field = (
                    ("beta", "Beta")
                    if key == "epigeneticAge"
                    else ("hazardRatio", "HR")
                )
                require_source_value(
                    record.get(effect_field[0]),
                    row.get(effect_field[1]),
                    f"{symbol} {key} effect",
                )
                require_source_value(
                    record.get("pValue"),
                    row.get("p"),
                    f"{symbol} {key} P value",
                )

        for record in evidence.get("transcriptomic", []):
            row = source_row(
                transcriptomic,
                record["sourceRow"],
                2,
                f"{symbol} transcriptomics",
            )
            for field, source_field in (
                ("sourceSymbol", "Gene.symbol"),
                ("slope", "Slope"),
                ("standardError", "SE"),
                ("pearsonCorrelation", "Pearson.corr"),
                ("pValue", "P.Value"),
                ("adjustedPValue", "P.Adjusted"),
            ):
                require_source_value(
                    record.get(field),
                    row.get(source_field),
                    f"{symbol} transcriptomics {field}",
                )

        for record in evidence.get("organAge", []):
            row = source_row(
                organage,
                record["sourceRow"],
                2,
                f"{symbol} OrganAge",
            )
            for field, source_field in (
                ("organ", "organ"),
                ("seqId", "seq_id"),
                ("sourceSymbol", "gene_symbol"),
                ("selectedModels", "selected_models"),
                ("meanNonzeroCoefficient", "mean_nonzero_coefficient"),
                ("sourceCommit", "source_commit"),
            ):
                require_source_value(
                    record.get(field),
                    row.get(source_field),
                    f"{symbol} OrganAge {field}",
                )
            if record.get("sourceCommit") != ORGANAGE_COMMIT:
                fail(f"{symbol}: unpinned OrganAge source")

    for protein_id, protein in entities["protein"].items():
        source = uniprot_rows.get(protein_id)
        if not source:
            fail(f"{protein_id}: missing from UniProtKB source file")
        require_source_value(
            protein.get("entryName"),
            source.get("Entry Name"),
            f"{protein_id} UniProtKB entry name",
        )
        require_source_value(
            protein.get("function"),
            uniprot_function(source.get("Function [CC]")),
            f"{protein_id} UniProtKB function",
        )

    for cpg_id, cpg in entities["cpg"].items():
        for item in cpg.get("sourceDescriptions", []):
            require_source_value(
                item.get("text"),
                epigenetic_titles.get(item.get("sourceSheet")),
                f"{cpg_id} source description",
            )

    for metabolite_id, metabolite in entities["metabolite"].items():
        source = chebi_rows.get(metabolite_id)
        if not source:
            fail(f"{metabolite_id}: missing from ChEBI source file")
        require_source_value(
            metabolite.get("name"),
            plain_text(source.get("name")),
            f"{metabolite_id} ChEBI name",
        )
        require_source_value(
            metabolite.get("description"),
            plain_text(source.get("definition")),
            f"{metabolite_id} ChEBI definition",
        )
        for record in metabolite.get("agingEvidence", []):
            row = source_row(
                metaboage,
                record["sourceRow"],
                2,
                f"{metabolite_id} MetaboAge",
            )
            for field, source_field in (
                ("metaboliteName", "Metabolite"),
                ("method", "Method"),
                ("value", "Mean value/beta value"),
                ("pmid", "Study  PMID reference"),
            ):
                expected = row.get(source_field)
                if field == "pmid":
                    pmid_match = re.search(r"\d+", normalized_text(expected))
                    expected = pmid_match.group(0) if pmid_match else None
                require_source_value(
                    record.get(field),
                    expected,
                    f"{metabolite_id} MetaboAge {field}",
                )

    for pathway_id, pathway in entities["pathway"].items():
        snapshot = reactome_snapshot.get(pathway_id)
        if not snapshot:
            fail(f"{pathway_id}: missing Reactome description snapshot")
        expected_description = (
            plain_text(snapshot.get("summary"))
            or plain_text(snapshot.get("name"))
        )
        require_source_value(
            pathway.get("description"),
            expected_description,
            f"{pathway_id} Reactome description",
        )
        for entity_type, members, source_file in (
            ("gene", pathway.get("genes", []), "NCBI2Reactome_PE_Pathway.txt"),
            (
                "protein",
                pathway.get("proteins", []),
                "UniProt2Reactome_PE_Pathway.txt",
            ),
            (
                "metabolite",
                pathway.get("metabolites", []),
                "ChEBI2Reactome_PE_Pathway.txt",
            ),
        ):
            for member in members:
                row = reactome_rows[source_file].get(member["sourceRow"])
                if not row or len(row) < 8:
                    fail(
                        f"{pathway_id} {entity_type}: invalid Reactome "
                        f"source row {member.get('sourceRow')}"
                    )
                expected_id = (
                    member["sourceIdentifier"]
                    if entity_type == "gene"
                    else member["id"]
                )
                source_id = row[0]
                if entity_type == "metabolite" and not source_id.startswith(
                    "CHEBI:"
                ):
                    source_id = f"CHEBI:{source_id}"
                require_source_value(
                    expected_id,
                    source_id,
                    f"{pathway_id} {entity_type} identifier",
                )
                require_source_value(
                    member.get("physicalEntityId"),
                    row[1],
                    f"{pathway_id} {entity_type} physical entity",
                )
                require_source_value(
                    pathway_id,
                    row[3],
                    f"{pathway_id} {entity_type} pathway",
                )
                require_source_value(
                    row[7],
                    "Homo sapiens",
                    f"{pathway_id} {entity_type} species",
                )
def load_entities(manifest: dict[str, Any]) -> dict[str, dict[str, dict[str, Any]]]:
    entities: dict[str, dict[str, dict[str, Any]]] = {}
    for entity_type, config in manifest["entities"].items():
        records: dict[str, dict[str, Any]] = {}
        for chunk in config["chunks"]:
            records.update(load_json(DATA_DIR / "chunks" / chunk))
        if len(records) != config["count"]:
            fail(f"{entity_type}: manifest count differs from chunks")
        if set(records) != set(config["idToChunk"]):
            fail(f"{entity_type}: id-to-chunk index differs from chunk records")
        entities[entity_type] = records
    return entities


def validate_gene_evidence(symbol: str, gene: dict[str, Any]) -> None:
    evidence = gene.get("evidence", {})
    genage = evidence.get("genAgeHuman")
    if genage:
        require_fields(
            genage,
            ("genAgeId", "evidenceBasis", "sourceUrl"),
            f"{symbol} GenAge",
        )
        require_url(genage["sourceUrl"], f"{symbol} GenAge")
    for row in evidence.get("longevityMap", []):
        require_fields(
            row,
            ("reportId", "association", "sourceUrl"),
            f"{symbol} LongevityMap",
        )
        require_url(row["sourceUrl"], f"{symbol} LongevityMap")
        if row.get("pubmedUrl"):
            require_url(row["pubmedUrl"], f"{symbol} LongevityMap publication")
    for key in ("epigeneticAge", "epigeneticMortality"):
        for row in evidence.get(key, []):
            require_fields(
                row,
                (
                    "cpg",
                    "endpoint",
                    "model",
                    "pValue",
                    "sourceSheet",
                    "sourceRow",
                    "sourceUrl",
                    "publicationUrl",
                ),
                f"{symbol} {key}",
            )
            require_url(row["sourceUrl"], f"{symbol} {key}")
            require_url(
                row["publicationUrl"],
                f"{symbol} {key} publication",
            )
    for row in evidence.get("transcriptomic", []):
        require_fields(
            row,
            (
                "organism",
                "cohort",
                "endpoint",
                "model",
                "slope",
                "adjustedPValue",
                "sourceSymbol",
                "sourceSheet",
                "sourceRow",
                "sourceUrl",
                "publicationUrl",
            ),
            f"{symbol} transcriptomics",
        )
        if row["organism"] != "Human":
            fail(f"{symbol}: non-human transcriptomic record retained")
        if any(key.lower().startswith("sourcemouse") for key in row):
            fail(f"{symbol}: cross-species mapping field retained")
        require_url(row["sourceUrl"], f"{symbol} transcriptomics")
        require_url(
            row["publicationUrl"],
            f"{symbol} transcriptomics publication",
        )
    for row in evidence.get("organAge", []):
        require_fields(
            row,
            (
                "organism",
                "organ",
                "targetName",
                "selectedModels",
                "sourceCommit",
                "sourceRow",
                "sourceUrl",
                "publicationUrl",
            ),
            f"{symbol} OrganAge",
        )
        if row["organism"] != "Human":
            fail(f"{symbol}: non-human OrganAge record retained")
        require_url(row["sourceUrl"], f"{symbol} OrganAge")
        require_url(
            row["publicationUrl"],
            f"{symbol} OrganAge publication",
        )


def validate() -> dict[str, int]:
    manifest = load_json(DATA_DIR / "manifest.json")
    collections = load_json(DATA_DIR / "collections.json")
    search_index = load_json(DATA_DIR / "search-index.json")
    quality = load_json(DATA_DIR / "quality-report.json")
    audit = load_json(DATA_DIR / "selection-audit.json")
    entities = load_entities(manifest)

    if set(entities) != EXPECTED_TYPES:
        fail(f"entity types are {set(entities)}, expected {EXPECTED_TYPES}")
    module_keys = {row["key"] for row in collections["modules"]}
    if module_keys != EXPECTED_MODULES:
        fail(f"module keys are {module_keys}, expected {EXPECTED_MODULES}")
    if any("reaction" in path.name.lower() for path in (DATA_DIR / "chunks").iterdir()):
        fail("reaction chunk exists in reduced release")
    if len(entities["gene"]) != 1000:
        fail("release must contain exactly 1,000 genes")
    if len(audit) != 1000:
        fail("selection audit must contain exactly 1,000 genes")
    if [row["releaseRank"] for row in audit] != list(range(1, 1001)):
        fail("selection audit ranks are not complete")

    index_keys = {(row["type"], row["id"]) for row in search_index}
    entity_keys = {
        (entity_type, entity_id)
        for entity_type, records in entities.items()
        for entity_id in records
    }
    if index_keys != entity_keys:
        fail("search index and entity records differ")

    for entity_type in (
        "gene",
        "cpg",
        "protein",
        "metabolite",
        "pathway",
    ):
        coverage = quality["sourceDescriptionCoverage"][entity_type]
        if coverage["withDescription"] != coverage["total"]:
            fail(f"{entity_type}: incomplete source description coverage")

    genes = entities["gene"]
    proteins = entities["protein"]
    cpgs = entities["cpg"]
    metabolites = entities["metabolite"]
    pathways = entities["pathway"]

    for symbol, gene in genes.items():
        require_fields(
            gene,
            (
                "symbol",
                "name",
                "description",
                "descriptionSource",
                "summarySource",
                "annotation",
                "evidenceLayers",
                "connections",
            ),
            symbol,
        )
        source = gene["summarySource"]
        require_fields(source, ("label", "url"), f"{symbol} description source")
        if source["label"] not in {
            "NCBI Gene summary",
            "HGNC approved gene description",
        }:
            fail(f"{symbol}: unsupported gene description source")
        require_url(source["url"], f"{symbol} description")
        annotation = gene["annotation"]
        require_fields(
            annotation,
            (
                "approvedName",
                "hgncId",
                "humanEntrezId",
                "chromosomeLocation",
                "hgncUrl",
            ),
            f"{symbol} annotation",
        )
        require_url(annotation["hgncUrl"], f"{symbol} HGNC")
        validate_gene_evidence(symbol, gene)
        connections = gene["connections"]
        for protein_id in connections.get("proteins", []):
            if protein_id not in proteins:
                fail(f"{symbol}: unresolved protein {protein_id}")
            if proteins[protein_id].get("geneSymbol") != symbol:
                fail(f"{symbol}: protein {protein_id} maps to another gene")
        for cpg_id in connections.get("cpgs", []):
            if cpg_id not in cpgs:
                fail(f"{symbol}: unresolved CpG {cpg_id}")
            mapped = {row["symbol"] for row in cpgs[cpg_id]["geneMappings"]}
            if symbol not in mapped:
                fail(f"{symbol}: CpG {cpg_id} lacks reciprocal mapping")
        for pathway in connections.get("pathways", []):
            pathway_id = pathway["id"]
            if pathway_id not in pathways:
                fail(f"{symbol}: unresolved pathway {pathway_id}")
            require_fields(
                pathway,
                (
                    "physicalEntityId",
                    "mappingScope",
                    "sourceFile",
                    "sourceRow",
                    "sourceUrl",
                    "pathwayUrl",
                ),
                f"{symbol} pathway {pathway_id}",
            )
            require_url(pathway["sourceUrl"], f"{symbol} pathway entity")
            members = {row["symbol"] for row in pathways[pathway_id]["genes"]}
            if symbol not in members:
                fail(f"{symbol}: pathway {pathway_id} lacks reciprocal gene")

    for cpg_id, cpg in cpgs.items():
        require_fields(
            cpg,
            (
                "name",
                "description",
                "descriptionSource",
                "sourceDescriptions",
                "geneMappings",
                "agingEvidence",
                "url",
            ),
            cpg_id,
        )
        require_url(cpg["url"], cpg_id)
        for item in cpg["sourceDescriptions"]:
            require_fields(
                item,
                (
                    "source",
                    "sourceSheet",
                    "text",
                    "sourceUrl",
                    "publicationUrl",
                ),
                f"{cpg_id} source description",
            )
            require_url(
                item["sourceUrl"],
                f"{cpg_id} source description table",
            )
            require_url(
                item["publicationUrl"],
                f"{cpg_id} source description publication",
            )
        for mapping in cpg["geneMappings"]:
            if mapping["symbol"] not in genes:
                fail(f"{cpg_id}: unresolved gene {mapping['symbol']}")
            require_url(mapping["sourceUrl"], f"{cpg_id} gene mapping")
        for row in cpg["agingEvidence"]:
            require_fields(
                row,
                (
                    "source",
                    "endpoint",
                    "model",
                    "effectType",
                    "effect",
                    "pValue",
                    "sourceSheet",
                    "sourceRow",
                    "sourceUrl",
                    "publicationUrl",
                ),
                f"{cpg_id} evidence",
            )
            require_url(row["sourceUrl"], f"{cpg_id} evidence")
            require_url(
                row["publicationUrl"],
                f"{cpg_id} evidence publication",
            )

    for protein_id, protein in proteins.items():
        require_fields(
            protein,
            (
                "name",
                "geneSymbol",
                "description",
                "descriptionSource",
                "url",
                "organAgeEvidence",
            ),
            protein_id,
        )
        require_url(protein["url"], protein_id)
        if protein["descriptionSource"] != "UniProtKB Function annotation":
            fail(f"{protein_id}: protein lacks a UniProtKB function annotation")
        if protein["geneSymbol"] not in genes:
            fail(f"{protein_id}: unresolved encoding gene")
        for pathway in protein.get("pathways", []):
            if pathway["id"] not in pathways:
                fail(f"{protein_id}: unresolved pathway {pathway['id']}")
            require_url(pathway["sourceUrl"], f"{protein_id} pathway entity")
            pathway_proteins = {
                row["id"] for row in pathways[pathway["id"]]["proteins"]
            }
            if protein_id not in pathway_proteins:
                fail(f"{protein_id}: pathway lacks reciprocal protein")
        for row in protein.get("organAgeEvidence", []):
            require_url(row["sourceUrl"], f"{protein_id} OrganAge")
            require_url(
                row.get("publicationUrl"),
                f"{protein_id} OrganAge publication",
            )

    for metabolite_id, metabolite in metabolites.items():
        require_fields(
            metabolite,
            (
                "name",
                "description",
                "descriptionSource",
                "url",
                "agingEvidence",
            ),
            metabolite_id,
        )
        require_url(metabolite["url"], metabolite_id)
        if metabolite["descriptionSource"] != "ChEBI definition":
            fail(f"{metabolite_id}: metabolite lacks a ChEBI definition")
        for pathway in metabolite.get("pathways", []):
            if pathway["id"] not in pathways:
                fail(f"{metabolite_id}: unresolved pathway {pathway['id']}")
            require_url(pathway["sourceUrl"], f"{metabolite_id} pathway entity")
            pathway_metabolites = {
                row["id"] for row in pathways[pathway["id"]]["metabolites"]
            }
            if metabolite_id not in pathway_metabolites:
                fail(f"{metabolite_id}: pathway lacks reciprocal metabolite")
        for row in metabolite.get("agingEvidence", []):
            require_fields(
                row,
                (
                    "source",
                    "metaboliteName",
                    "sourceFile",
                    "sourceSheet",
                    "sourceRow",
                    "sourceUrl",
                ),
                f"{metabolite_id} MetaboAge",
            )
            require_url(row["sourceUrl"], f"{metabolite_id} MetaboAge")
            if row.get("pmid"):
                require_url(
                    row.get("publicationUrl"),
                    f"{metabolite_id} MetaboAge publication",
                )

    for pathway_id, pathway in pathways.items():
        require_fields(
            pathway,
            (
                "name",
                "description",
                "descriptionSource",
                "source",
                "species",
                "url",
                "recordUrl",
                "genes",
            ),
            pathway_id,
        )
        if pathway["source"] != "Reactome" or pathway["species"] != "Homo sapiens":
            fail(f"{pathway_id}: invalid pathway source or species")
        if pathway["descriptionSource"] not in {
            "Reactome summation",
            "Reactome pathway name",
        }:
            fail(f"{pathway_id}: invalid description source")
        require_url(pathway["url"], pathway_id)
        require_url(pathway["recordUrl"], f"{pathway_id} record")
        for member in pathway["genes"]:
            if member["symbol"] not in genes:
                fail(f"{pathway_id}: unresolved gene {member['symbol']}")
            if member.get("mappingScope") != "Direct Reactome physical-entity mapping":
                fail(f"{pathway_id}: indirect gene mapping retained")
            require_url(member["sourceUrl"], f"{pathway_id} gene")
        for member in pathway.get("proteins", []):
            if member["id"] not in proteins:
                fail(f"{pathway_id}: unresolved protein {member['id']}")
            require_url(member["sourceUrl"], f"{pathway_id} protein")
        for member in pathway.get("metabolites", []):
            if member["id"] not in metabolites:
                fail(f"{pathway_id}: unresolved metabolite {member['id']}")
            require_url(member["sourceUrl"], f"{pathway_id} metabolite")

    for record in walk(entities):
        if "reaction" in {key.lower() for key in record}:
            fail("reaction field retained in published entity")
        for value in record.values():
            if isinstance(value, str) and any(
                phrase.lower() in value.lower()
                for phrase in FORBIDDEN_GENERATED_PHRASES
            ):
                fail(f"generated placeholder phrase retained: {value[:100]}")

    validate_source_replay(entities)

    return {
        "genes": len(genes),
        "cpgs": len(cpgs),
        "proteins": len(proteins),
        "metabolites": len(metabolites),
        "pathways": len(pathways),
        "searchRecords": len(search_index),
    }


if __name__ == "__main__":
    print(json.dumps(validate(), indent=2))
