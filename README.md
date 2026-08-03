# Human Aging Atlas

[Open the Human Aging Atlas](https://moqri-lab.org/Human-Aging-Atlas/)

The Human Aging Atlas is a source-traceable research resource developed by the
[Moqri Lab](https://moqri.bwh.harvard.edu/). It brings together human aging and
longevity evidence across several biological layers and links every displayed
record to its originating dataset or reference resource.

The current community preview **(V3)** is a static release designed for rapid,
reproducible browsing. It does not require a database and does not call external
APIs while a user is browsing the site.

## Evidence included

- **Genomics:** GenAge and LongevityMap
- **Epigenomics:** chronological-age and mortality-associated CpG records
- **Transcriptomics:** human multi-tissue age-associated expression records
- **Proteomics:** OrganAge evidence with UniProtKB protein identity
- **Metabolomics:** MetaboAge DB evidence with ChEBI identity
- **Pathways:** direct human physical-entity mappings from Reactome

Gene, protein, and metabolite pathway links require direct human Reactome
physical-entity mappings. Shared pathway membership is not converted into a
gene-metabolite relationship, and pathway context is not counted as an
independent aging association.

## Source descriptions

- Gene descriptions come from NCBI Gene, with the HGNC approved name used when
  an NCBI summary is unavailable.
- CpG descriptions use the exact titles of the corresponding source tables.
- Protein descriptions come from UniProtKB Function annotations.
- Metabolite descriptions come from ChEBI.
- Pathway descriptions come from Reactome ContentService summations or the
  exact Reactome pathway name.

Records without the required source-supported identity or description are not
published. The Atlas does not use generated scientific prose for these fields.

## Data and provenance

The deployed website is self-contained in `index.html`, `app.js`, `styles.css`,
`assets/`, and `data/`. Generated data are divided into small JSON chunks for
fast static loading.

Important release records include:

- `data/manifest.json`: release rules, entity index, and source checksums
- `data/selection-audit.json`: deterministic gene-selection audit
- `data/quality-report.json`: field and source-description coverage
- `data/search-index.json`: searchable Atlas index

Raw-source locations and loaders are defined in `build/source_loaders.py`.
Source files that cannot be redistributed are not included in this repository;
their checksums and source identities remain recorded in the release manifest.

## Rebuild and validate

After the required source files are available, run:

```bash
python3 build/fetch_source_descriptions.py
python3 build/build_data.py
python3 build/validate_data.py
```

Validation replays retained GenAge, LongevityMap, epigenomic, transcriptomic,
OrganAge, and MetaboAge values against their local source rows. It also checks
UniProtKB functions, ChEBI definitions, Reactome descriptions, and direct
Reactome mapping rows. A source-row mismatch stops the build.

## Run locally

Double-click `Open Human Aging Atlas.command`, or serve the repository with a
local static HTTP server. GitHub Actions automatically publishes the `main`
branch to GitHub Pages.

## Project and attribution

Developed by the [Moqri Lab](https://moqri.bwh.harvard.edu/).

Current community preview implementation and data integration:
**Rey Zafarnejad, PhD, MS**.

## Use

This preview is intended for research and scientific communication. It does not
provide clinical recommendations or establish causality.
