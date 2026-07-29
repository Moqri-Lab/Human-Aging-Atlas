# Human Aging Atlas V3

Version 3 is a source-traceable static preview of the Human Aging Atlas. It is
kept separate from earlier project deliverables.

The current static release contains 1,000 human genes and six searchable
evidence areas:

- Genomics: GenAge and LongevityMap
- Epigenomics: cAge and bAge
- Transcriptomics: human tAge records
- Proteomics: OrganAge evidence with UniProtKB protein identity
- Metabolomics: MetaboAge DB evidence with ChEBI identity
- Pathways: direct human Reactome physical-entity mappings

Gene, protein, and metabolite pathway links come directly from the
corresponding Reactome physical-entity mapping files. Shared pathway membership
is not converted into a gene-metabolite relationship.

Gene descriptions come from NCBI Gene and fall back to the exact HGNC approved
name when NCBI does not provide a summary. CpG descriptions are the exact
titles of the applicable source workbook tables. Protein descriptions come
from UniProtKB Function annotations; proteins without that annotation are not
published. Metabolite descriptions come from ChEBI; metabolites without a
ChEBI definition are not published. Pathway descriptions come from Reactome
ContentService summations or the exact Reactome pathway name.

## Rebuild and validate

Use the workspace Python runtime:

```bash
python3 build/fetch_source_descriptions.py
python3 build/build_data.py
python3 build/validate_data.py
```

The first command refreshes the local NCBI Gene and Reactome description
snapshots. The website itself does not call external APIs at runtime.

The build writes:

- `data/selection-audit.json`: deterministic 1,000-gene selection
- `data/quality-report.json`: field and description coverage
- `data/manifest.json`: release rules, source checksums, and chunk index

Validation replays retained GenAge, LongevityMap, cAge/bAge, tAge, OrganAge,
and MetaboAge values against their local source rows. It also checks UniProtKB
functions, ChEBI definitions, Reactome descriptions, and direct Reactome
mapping rows. A source-row mismatch stops the build.

## Website

The GitHub Pages deployment serves the static files in this repository. The
website does not require a database or call external APIs at runtime.

## Run locally

Double-click `Open Human Aging Atlas.command`, or serve this folder with a
local static HTTP server. The production backend can later generate the same
records from normalized Parquet tables and PostgreSQL.
