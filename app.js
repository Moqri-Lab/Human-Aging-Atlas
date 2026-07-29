const app = document.getElementById("app");
const DATA_VERSION = "14";

const ENTITY_LABELS = {
  gene: "Gene",
  cpg: "CpG site",
  protein: "Protein",
  metabolite: "Metabolite",
  pathway: "Pathway",
};

const ENTITY_MODULES = {
  gene: "genomics",
  cpg: "epigenomics",
  protein: "proteomics",
  metabolite: "metabolomics",
  pathway: "pathways",
};

const MODULE_TONES = {
  genomics: "gold",
  epigenomics: "blue",
  transcriptomics: "teal",
  proteomics: "slate",
  metabolomics: "green",
  pathways: "navy",
};

const state = {
  manifest: null,
  collections: null,
  searchIndex: [],
  chunkCache: new Map(),
  atlas: { query: "", sort: "gene", direction: "asc", page: 1 },
  module: { key: "", query: "", sort: "id", direction: "asc", page: 1 },
};

let tableSequence = 0;

function escapeHtml(value) {
  return String(value ?? "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}

function safeUrl(value) {
  if (!value) return "";
  try {
    const url = new URL(value, window.location.href);
    return ["http:", "https:"].includes(url.protocol) ? url.href : "";
  } catch {
    return "";
  }
}

function externalLink(url, label, className = "") {
  const safe = safeUrl(url);
  return safe
    ? `<a class="${className}" href="${escapeHtml(safe)}" target="_blank" rel="noreferrer">${escapeHtml(label)}</a>`
    : escapeHtml(label);
}

function sourceLink(url, label = "Verify") {
  const safe = safeUrl(url);
  return safe
    ? `<a class="source-link" href="${escapeHtml(safe)}" target="_blank" rel="noreferrer">${escapeHtml(label)} ↗</a>`
    : "";
}

function sourceTextBlock(text, url, sourceLabel) {
  const value = String(text ?? "").trim();
  if (!value) return "";
  const source = sourceLink(url, sourceLabel);
  if (value.length <= 850) {
    return `<p>${escapeHtml(value)}</p>${source}`;
  }
  const candidate = value.slice(0, 700);
  const sentenceEnd = Math.max(
    candidate.lastIndexOf(". "),
    candidate.lastIndexOf("? "),
    candidate.lastIndexOf("! "),
  );
  const excerpt =
    sentenceEnd >= 350
      ? value.slice(0, sentenceEnd + 1)
      : `${candidate.trimEnd()}…`;
  return `
    <p>${escapeHtml(excerpt)}</p>
    <details class="source-text-more">
      <summary>Read complete source text</summary>
      <p>${escapeHtml(value)}</p>
    </details>
    ${source}`;
}

function sourcedValue(value, url) {
  const safe = safeUrl(url);
  return safe
    ? `<a class="sourced-value" href="${escapeHtml(safe)}" target="_blank" rel="noreferrer">${value}</a>`
    : value;
}

function normalize(value) {
  return String(value ?? "")
    .toLowerCase()
    .normalize("NFKD")
    .replace(/[^\w]+/g, " ")
    .trim();
}

function formatValue(value) {
  if (value === null || value === undefined || value === "") return "—";
  if (typeof value === "object" && value.display) {
    return escapeHtml(value.display);
  }
  return escapeHtml(value);
}

function scientific(value) {
  if (value === null || value === undefined || value === "") return "—";
  if (typeof value === "object" && value.display) {
    return escapeHtml(value.display);
  }
  const number = Number(typeof value === "object" ? value.value : value);
  if (!Number.isFinite(number)) return escapeHtml(value);
  if (number === 0) return "0";
  return Math.abs(number) < 0.001
    ? number.toExponential(2)
    : Number(number.toPrecision(4)).toString();
}

function naturalCompare(left, right) {
  return String(left).localeCompare(String(right), undefined, {
    numeric: true,
    sensitivity: "base",
  });
}

function hashPath(type, id, focus = "") {
  const suffix = focus ? `/${encodeURIComponent(focus)}` : "";
  return `#/${type}/${encodeURIComponent(id)}${suffix}`;
}

function entityLink(type, id, label = id, className = "", focus = "") {
  return `<a class="${className}" href="${hashPath(type, id, focus)}">${escapeHtml(label)}</a>`;
}

function indexedEntityLink(type, id, label = id, className = "") {
  return state.manifest.entities[type]?.idToChunk?.[id]
    ? entityLink(type, id, label, className)
    : escapeHtml(label);
}

async function getJson(path) {
  const separator = path.includes("?") ? "&" : "?";
  const response = await fetch(`${path}${separator}v=${DATA_VERSION}`);
  if (!response.ok) throw new Error(`Unable to load ${path}`);
  return response.json();
}

async function loadEntity(type, id) {
  const config = state.manifest.entities[type];
  const chunk = config?.idToChunk?.[id];
  if (!chunk) return null;
  const cacheKey = `${type}:${chunk}`;
  if (!state.chunkCache.has(cacheKey)) {
    state.chunkCache.set(cacheKey, await getJson(`data/chunks/${chunk}`));
  }
  return state.chunkCache.get(cacheKey)[id] ?? null;
}

function parseRoute() {
  const path = window.location.hash.replace(/^#\/?/, "") || "atlas";
  const [section, rawId = "", rawFocus = ""] = path.split("/");
  return {
    section,
    id: decodeURIComponent(rawId),
    focus: decodeURIComponent(rawFocus),
  };
}

function moduleDefinition(key) {
  return state.collections.modules.find((item) => item.key === key);
}

function navLink(route, label, active) {
  return `<a href="#/${route}" class="top-nav-link ${active ? "active" : ""}">${escapeHtml(label)}</a>`;
}

function evidenceBadges(layers = [], compact = false) {
  if (!layers.length) return `<span class="muted">No direct aging evidence</span>`;
  return `<div class="evidence-badges ${compact ? "compact" : ""}">
    ${layers
      .map((item) => {
        const module = moduleDefinition(item.key) ?? {
          label: item.key,
          key: item.key,
        };
        const sources =
          !compact && item.sources?.length
            ? ` (${item.sources.join(" + ")})`
            : "";
        return `<a href="#/module/${encodeURIComponent(item.key)}" class="badge ${MODULE_TONES[item.key] ?? "slate"}">${escapeHtml(module.label + sources)}</a>`;
      })
      .join("")}
  </div>`;
}

function compactSearch(activeSection) {
  if (activeSection === "atlas" || ENTITY_LABELS[activeSection]) return "";
  return `
    <form class="compact-search" id="compact-search">
      <input id="compact-search-input" type="search" aria-label="Search atlas" placeholder="Search the atlas" />
      <button type="submit" class="icon-search" aria-label="Open first match" title="Open first match">⌕</button>
    </form>`;
}

function entityPageSearch(activeSection) {
  if (!ENTITY_LABELS[activeSection]) return "";
  return `
    <form class="entity-page-search" id="entity-page-search">
      <label for="entity-page-search-input">Search the atlas</label>
      <div>
        <input id="entity-page-search-input" type="search" autocomplete="off" placeholder="Gene, protein, metabolite, CpG site, or pathway" />
        <button type="submit">Open first match</button>
      </div>
    </form>`;
}

function contextNavigation(activeSection) {
  if (activeSection !== "module" && !ENTITY_LABELS[activeSection]) return "";
  return `
    <nav class="context-nav" aria-label="Page navigation">
      <button type="button" class="context-back" id="context-back" aria-label="Go back">
        <span aria-hidden="true">←</span><span>Back</span>
      </button>
      <a href="#/atlas" class="context-home">Atlas home</a>
    </nav>`;
}

function sharedRail(activeSection, activeModule = "") {
  const links = state.collections.modules
    .map(
      (module) => `
        <a href="#/module/${module.key}" class="rail-link ${activeModule === module.key ? "active" : ""}">
          <strong>${escapeHtml(module.label)}</strong>
        </a>`,
    )
    .join("");
  return `
    <aside class="rail" aria-label="Atlas navigation">
      <a href="#/atlas" class="rail-home ${activeSection === "atlas" ? "active" : ""}">
        <strong>Atlas</strong>
      </a>
      <div class="rail-heading">Evidence</div>
      <nav class="rail-group">${links}</nav>
    </aside>`;
}

function shell(content, activeSection = "atlas", activeModule = "") {
  const atlasActive = ["atlas", "module", ...Object.keys(ENTITY_LABELS)].includes(
    activeSection,
  );
  app.innerHTML = `
    <header class="site-header">
      <a href="#/atlas" class="brand">
        <span class="brand-mark"><img src="assets/atlas-logo.png" alt="" /></span>
        <span>Human Aging Atlas</span>
      </a>
      <nav class="top-nav" aria-label="Primary navigation">
        ${navLink("atlas", "Atlas", atlasActive)}
        ${navLink("methods", "Methods", activeSection === "methods")}
        ${navLink("sources", "Sources", activeSection === "sources")}
      </nav>
      ${compactSearch(activeSection)}
    </header>
    <div class="page-grid">
      ${sharedRail(activeSection, activeModule)}
      <main class="main-content">
        ${contextNavigation(activeSection)}
        ${entityPageSearch(activeSection)}
        ${content}
      </main>
    </div>
    <footer class="site-footer">
      <span>Human Aging Atlas</span>
      <span>Source-traceable human aging evidence.</span>
    </footer>`;
  bindNavigation();
  bindSearchForms();
  bindRevealers();
}

function bindNavigation() {
  document.getElementById("context-back")?.addEventListener("click", () => {
    if (window.history.length > 1) window.history.back();
    else window.location.hash = "#/atlas";
  });
}

function bindSearchForms() {
  const forms = [
    ["compact-search", "compact-search-input"],
    ["entity-page-search", "entity-page-search-input"],
  ];
  forms.forEach(([formId, inputId]) => {
    document.getElementById(formId)?.addEventListener("submit", (event) => {
      event.preventDefault();
      const query = document.getElementById(inputId).value.trim();
      const result = searchRecords(query)[0];
      if (result) window.location.hash = hashPath(result.type, result.id);
    });
  });
}

function pageHeader(eyebrow, title, description = "") {
  return `
    <header class="page-header">
      <div class="eyebrow">${escapeHtml(eyebrow)}</div>
      <h1>${escapeHtml(title)}</h1>
      ${description ? `<p>${escapeHtml(description)}</p>` : ""}
    </header>`;
}

function searchRecords(query, moduleKey = "") {
  const needle = normalize(query);
  const module = moduleKey ? moduleDefinition(moduleKey) : null;
  return state.searchIndex
    .filter((record) => {
      if (moduleKey && !record.modules?.includes(moduleKey)) return false;
      if (module?.primaryType && record.type !== module.primaryType) return false;
      if (!needle) return true;
      const haystack = normalize(
        [
          record.id,
          record.name,
          record.description,
          ...(record.aliases ?? []),
        ].join(" "),
      );
      return haystack.includes(needle);
    })
    .sort((left, right) => {
      if (!needle) return naturalCompare(left.id, right.id);
      const leftId = normalize(left.id);
      const rightId = normalize(right.id);
      const leftRank = leftId === needle ? 0 : leftId.startsWith(needle) ? 1 : 2;
      const rightRank = rightId === needle ? 0 : rightId.startsWith(needle) ? 1 : 2;
      return leftRank - rightRank || naturalCompare(left.id, right.id);
    });
}

function tableHtml(headers, rows, options = {}) {
  if (!rows.length) return "";
  const limit = options.limit ?? 20;
  const tableId = `table-${++tableSequence}`;
  return `
    <div class="table-wrap">
      <table class="atlas-table columns-${headers.length}">
        <thead><tr>${headers.map((header) => `<th>${header}</th>`).join("")}</tr></thead>
        <tbody>
          ${rows
            .map(
              (row, index) => `
                <tr class="${index >= limit ? "reveal-row is-hidden" : ""}" data-table="${tableId}">
                  ${row.map((cell) => `<td>${cell ?? ""}</td>`).join("")}
                </tr>`,
            )
            .join("")}
        </tbody>
      </table>
    </div>
    ${
      rows.length > limit
        ? `<button class="show-more table-revealer" data-table="${tableId}">Show more</button>`
        : ""
    }`;
}

function bindRevealers() {
  document.querySelectorAll(".table-revealer").forEach((button) => {
    button.addEventListener("click", () => {
      document
        .querySelectorAll(`[data-table="${button.dataset.table}"]`)
        .forEach((row) => row.classList.remove("is-hidden"));
      button.remove();
    });
  });
}

function detailDisclosure(title, content, open = false, id = "") {
  if (!content) return "";
  return `
    <details class="detail-disclosure" ${open ? "open" : ""} ${id ? `id="${id}"` : ""}>
      <summary>${escapeHtml(title)}</summary>
      <div class="detail-disclosure-body">${content}</div>
    </details>`;
}

function identityHeader(typeLabel, id, name, url, sourceLabel) {
  return `
    <header class="entity-header">
      <div>
        <div class="eyebrow">${escapeHtml(typeLabel)}</div>
        <h1>${escapeHtml(id)}</h1>
        <p>${escapeHtml(name)}</p>
      </div>
      <div class="entity-source">
        <span class="eyebrow">Identifier source</span>
        ${sourceLink(url, sourceLabel)}
      </div>
    </header>`;
}

function sourceLinksForGene(record, layer) {
  const urls = record.details?.evidenceSourceUrls ?? {};
  if (layer === "transcriptomics") {
    const source = state.collections.sources.find((item) => item.name === "tAge");
    return [
      sourceLink(source?.publicationUrl, "Nature study"),
      sourceLink(urls.tAge ?? source?.url, "Supplementary Table 2 (XLSX)"),
    ]
      .filter(Boolean)
      .join(" ");
  }
  const sourceNames =
    layer === "genomics"
      ? ["GenAge", "LongevityMap"]
      : [];
  return sourceNames
    .filter((name) => urls[name])
    .map((name) => sourceLink(urls[name], name))
    .join(" ");
}

function renderAtlasGeneTable() {
  const host = document.getElementById("atlas-results");
  if (!host) return;
  const genes = searchRecords(state.atlas.query)
    .filter((record) => record.type === "gene")
    .sort((left, right) => {
      const comparison =
        state.atlas.sort === "evidence"
          ? (left.evidenceLayers?.length ?? 0) -
              (right.evidenceLayers?.length ?? 0) ||
            naturalCompare(left.id, right.id)
          : naturalCompare(left.id, right.id);
      return state.atlas.direction === "asc" ? comparison : -comparison;
    });
  const shown = genes.slice(0, state.atlas.page * 25);
  host.innerHTML = `
    <section class="gene-index">
      <div class="section-heading">
        <div><div class="eyebrow">Gene atlas</div><h2>Human aging genes</h2></div>
      </div>
      ${
        shown.length
          ? `<div class="table-wrap">
              <table class="atlas-table gene-atlas-table">
                <thead>
                  <tr>
                    <th><button class="sort-button" data-atlas-sort="gene">Human gene${state.atlas.sort === "gene" ? (state.atlas.direction === "asc" ? " ↑" : " ↓") : ""}</button></th>
                    <th>Approved name</th>
                    <th>Human locus</th>
                    <th><button class="sort-button" data-atlas-sort="evidence">Evidence layers${state.atlas.sort === "evidence" ? (state.atlas.direction === "asc" ? " ↑" : " ↓") : ""}</button></th>
                    <th>Description source</th>
                  </tr>
                </thead>
                <tbody>
                  ${shown
                    .map(
                      (record) => `
                        <tr>
                          <td>${entityLink("gene", record.id, record.id, "table-link")}</td>
                          <td>${escapeHtml(record.name)}</td>
                          <td>${formatValue(record.details?.location)}</td>
                          <td>${evidenceBadges(record.evidenceLayers)}</td>
                          <td>${sourceLink(record.details?.sourceUrl, "Source")}</td>
                        </tr>`,
                    )
                    .join("")}
                </tbody>
              </table>
            </div>
            ${shown.length < genes.length ? `<button class="show-more" id="atlas-show-more">Show more</button>` : ""}`
          : `<div class="empty-state">No indexed gene matches this search.</div>`
      }
    </section>`;
  host.querySelectorAll("[data-atlas-sort]").forEach((button) => {
    button.addEventListener("click", () => {
      const key = button.dataset.atlasSort;
      if (state.atlas.sort === key) {
        state.atlas.direction = state.atlas.direction === "asc" ? "desc" : "asc";
      } else {
        state.atlas.sort = key;
        state.atlas.direction = key === "evidence" ? "desc" : "asc";
      }
      renderAtlasGeneTable();
    });
  });
  document.getElementById("atlas-show-more")?.addEventListener("click", () => {
    state.atlas.page += 1;
    renderAtlasGeneTable();
  });
}

function renderAtlas() {
  state.atlas = { query: "", sort: "gene", direction: "asc", page: 1 };
  shell(
    `
      ${pageHeader(
        "Human source reference",
        "Human Aging Atlas",
        "Search human aging evidence and open each value at its source.",
      )}
      ${connectionsMap()}
      <section class="atlas-search-panel">
        <form id="atlas-search-form" class="hero-search">
          <label for="atlas-search-input">Search genes</label>
          <div class="hero-search-row">
            <input id="atlas-search-input" type="search" autocomplete="off" placeholder="TP53 or tumor protein p53" />
            <button type="submit">Open first gene</button>
          </div>
        </form>
      </section>
      <div id="atlas-results"></div>
    `,
    "atlas",
  );
  const form = document.getElementById("atlas-search-form");
  const input = document.getElementById("atlas-search-input");
  input.addEventListener("input", () => {
    state.atlas.query = input.value;
    state.atlas.page = 1;
    renderAtlasGeneTable();
  });
  form.addEventListener("submit", (event) => {
    event.preventDefault();
    const result = searchRecords(input.value).find(
      (record) => record.type === "gene",
    );
    if (result) window.location.hash = hashPath("gene", result.id);
  });
  renderAtlasGeneTable();
}

function moduleTableDefinition(moduleKey) {
  const sourceCell = (record) =>
    sourceLink(record.details?.sourceUrl, "Source");
  return {
    genomics: {
      headers: ["Human gene", "Approved name", "Location", "Evidence", "Verify"],
      row: (record) => [
        entityLink("gene", record.id, record.id, "table-link", "genomics"),
        escapeHtml(record.name),
        formatValue(record.details?.location),
        evidenceBadges(
          record.evidenceLayers.filter((item) => item.key === "genomics"),
        ),
        sourceLinksForGene(record, "genomics"),
      ],
    },
    epigenomics: {
      headers: ["CpG site", "Mapped gene", "Endpoint", "Verify"],
      row: (record) => [
        entityLink("cpg", record.id, record.id, "table-link"),
        (record.details?.geneSymbols ?? [])
          .map((symbol) => indexedEntityLink("gene", symbol, symbol))
          .join(", "),
        escapeHtml((record.details?.endpoints ?? []).join("; ")),
        sourceCell(record),
      ],
    },
    transcriptomics: {
      headers: ["Human gene", "Approved name", "Location", "Sources"],
      row: (record) => [
        entityLink("gene", record.id, record.id, "table-link", "transcriptomics"),
        escapeHtml(record.name),
        formatValue(record.details?.location),
        sourceLinksForGene(record, "transcriptomics"),
      ],
    },
    proteomics: {
      headers: ["UniProt accession", "Protein name", "Encoding gene", "Verify"],
      row: (record) => [
        entityLink("protein", record.id, record.id, "table-link"),
        escapeHtml(record.name),
        indexedEntityLink("gene", record.details?.geneSymbol, record.details?.geneSymbol),
        sourceCell(record),
      ],
    },
    metabolomics: {
      headers: ["ChEBI ID", "Metabolite name", "Verify"],
      row: (record) => [
        entityLink("metabolite", record.id, record.id, "table-link"),
        escapeHtml(record.name),
        sourceCell(record),
      ],
    },
    pathways: {
      headers: ["Reactome ID", "Pathway name", "Directly mapped entity types", "Verify"],
      row: (record) => [
        entityLink("pathway", record.id, record.id, "table-link"),
        escapeHtml(record.name),
        escapeHtml((record.details?.memberTypes ?? []).join(", ")),
        sourceCell(record),
      ],
    },
  }[moduleKey];
}

function renderModuleResults() {
  const host = document.getElementById("module-results");
  if (!host) return;
  const definition = moduleTableDefinition(state.module.key);
  const records = searchRecords(state.module.query, state.module.key).sort(
    (left, right) => {
      const comparison =
        state.module.sort === "name"
          ? naturalCompare(left.name, right.name)
          : naturalCompare(left.id, right.id);
      return state.module.direction === "asc" ? comparison : -comparison;
    },
  );
  const shown = records.slice(0, state.module.page * 25);
  host.innerHTML = shown.length
    ? `${tableHtml(definition.headers, shown.map(definition.row), { limit: shown.length })}
       ${shown.length < records.length ? `<button class="show-more" id="module-show-more">Show more</button>` : ""}`
    : `<div class="empty-state">No matching record was found.</div>`;
  bindRevealers();
  document.getElementById("module-show-more")?.addEventListener("click", () => {
    state.module.page += 1;
    renderModuleResults();
  });
}

function renderModule(key) {
  const module = moduleDefinition(key);
  if (!module) {
    renderNotFound("module", key);
    return;
  }
  state.module = { key, query: "", sort: "id", direction: "asc", page: 1 };
  shell(
    `
      ${pageHeader("Evidence", module.label, module.description)}
      <section class="module-search-panel">
        <label for="module-search-input">Search ${escapeHtml(module.label.toLowerCase())}</label>
        <input id="module-search-input" type="search" autocomplete="off" placeholder="Identifier or source name" />
      </section>
      <div id="module-results"></div>
    `,
    "module",
    key,
  );
  document.getElementById("module-search-input").addEventListener("input", (event) => {
    state.module.query = event.target.value;
    state.module.page = 1;
    renderModuleResults();
  });
  renderModuleResults();
}

function renderGene(record, focus = "") {
  const evidence = record.evidence ?? {};
  const connections = record.connections ?? {};
  const genAge = evidence.genAgeHuman
    ? [[
        escapeHtml((evidence.genAgeHuman.evidenceBasis ?? []).join(", ")),
        sourcedValue(formatValue(evidence.genAgeHuman.genAgeId), evidence.genAgeHuman.sourceUrl),
        sourceLink(evidence.genAgeHuman.sourceUrl, "GenAge"),
      ]]
    : [];
  const longevity = (evidence.longevityMap ?? []).map((item) => [
    formatValue(item.population),
    formatValue(item.variants),
    item.pubmedId
      ? externalLink(item.pubmedUrl, `PMID ${item.pubmedId}`)
      : "—",
    sourceLink(item.sourceUrl, `LongevityMap ${item.reportId}`),
  ]);
  const epigenomics = [
    ...(evidence.epigeneticAge ?? []),
    ...(evidence.epigeneticMortality ?? []),
  ].map((item) => [
    entityLink("cpg", item.cpg, item.cpg, "table-link"),
    escapeHtml(item.endpoint),
    item.beta !== undefined ? "Beta" : "Hazard ratio",
    sourcedValue(
      item.beta !== undefined
        ? scientific(item.beta)
        : scientific(item.hazardRatio),
      item.sourceUrl,
    ),
    sourcedValue(scientific(item.pValue), item.sourceUrl),
    sourceLink(item.sourceUrl, `${item.sourceSheet} source table`),
    sourceLink(item.publicationUrl, "Publication"),
  ]);
  const transcriptomics = (evidence.transcriptomic ?? []).map((item) => [
    sourcedValue(scientific(item.slope), item.sourceUrl),
    sourcedValue(scientific(item.standardError), item.sourceUrl),
    sourcedValue(scientific(item.pearsonCorrelation), item.sourceUrl),
    sourcedValue(scientific(item.pValue), item.sourceUrl),
    sourcedValue(scientific(item.adjustedPValue), item.sourceUrl),
    escapeHtml(item.direction),
    sourceLink(item.sourceUrl, "Supplementary Table 2 (XLSX)"),
    sourceLink(item.publicationUrl, "Nature study"),
  ]);
  const organAge = (evidence.organAge ?? []).map((item) => [
    escapeHtml(item.organ),
    escapeHtml(item.targetName),
    sourcedValue(scientific(item.meanNonzeroCoefficient), item.sourceUrl),
    escapeHtml(item.coefficientDirection),
    sourcedValue(formatValue(item.selectedModels), item.sourceUrl),
    sourceLink(item.sourceUrl, "Published models"),
    sourceLink(item.publicationUrl, "Publication"),
  ]);
  const proteinRows = (connections.proteins ?? []).map((accession) => [
    entityLink("protein", accession, accession, "table-link"),
    sourceLink(
      `https://www.uniprot.org/uniprotkb/${accession}/entry`,
      "UniProtKB",
    ),
  ]);
  const pathwayRows = (connections.pathways ?? []).map((item) => [
    entityLink("pathway", item.id, item.id, "table-link"),
    escapeHtml(item.name),
    escapeHtml(item.physicalEntity ?? ""),
    sourceLink(item.sourceUrl, "Reactome entity"),
    sourceLink(item.pathwayUrl, "Pathway"),
  ]);
  const cpgRows = (connections.cpgs ?? []).map((id) => [
    entityLink("cpg", id, id, "table-link"),
    sourceLink(
      "https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs13073-023-01161-y/MediaObjects/13073_2023_1161_MOESM4_ESM.xlsx",
      "cAge/bAge source table",
    ),
  ]);
  const summarySource = record.summarySource ?? {};
  const annotation = record.annotation ?? {};
  const idRows = [
    ["HGNC ID", sourcedValue(formatValue(annotation.hgncId), annotation.hgncUrl)],
    ["NCBI Gene ID", sourcedValue(formatValue(annotation.humanEntrezId), annotation.ncbiUrl)],
    [
      "Ensembl gene ID",
      annotation.ensemblGeneId
        ? sourcedValue(
            escapeHtml(annotation.ensemblGeneId),
            `https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=${encodeURIComponent(annotation.ensemblGeneId)}`,
          )
        : "—",
    ],
    ["Human locus", sourcedValue(formatValue(annotation.chromosomeLocation), annotation.hgncUrl)],
  ];
  const activeModule =
    focus && moduleDefinition(focus) ? focus : ENTITY_MODULES.gene;
  shell(
    `
      ${identityHeader("Human gene", record.symbol, record.name, annotation.hgncUrl, "HGNC")}
      <section class="gene-function">
        <div class="eyebrow">Source description</div>
        ${sourceTextBlock(record.summary, summarySource.url, summarySource.label)}
      </section>
      <section class="detail-stack">
        ${detailDisclosure("Identifiers", tableHtml(["Identifier type", "Value"], idRows))}
        ${detailDisclosure("Genomics: GenAge", tableHtml(["Evidence basis", "GenAge ID", "Verify"], genAge), false, "genomics")}
        ${detailDisclosure("Genomics: significant LongevityMap associations", tableHtml(["Population", "Variant", "Publication", "Verify"], longevity), false, "longevity")}
        ${detailDisclosure("Epigenomics", tableHtml(["CpG site", "Endpoint", "Effect measure", "Estimate", "P value", "Source table", "Publication"], epigenomics), false, "epigenomics")}
        ${detailDisclosure("Transcriptomics: human multi-tissue chronological-age model", tableHtml(["Slope", "Standard error", "Pearson correlation", "P value", "Adjusted P", "Direction", "Source data", "Study"], transcriptomics), false, "transcriptomics")}
        ${detailDisclosure("Proteomics: OrganAge", tableHtml(["Organ", "Target", "Mean nonzero coefficient", "Direction", "Selected models", "Published models", "Publication"], organAge), false, "proteomics")}
        ${detailDisclosure("Encoded proteins", tableHtml(["UniProt accession", "Verify"], proteinRows))}
        ${detailDisclosure("Direct Reactome pathways", tableHtml(["Reactome ID", "Pathway", "Mapped entity", "Verify entity", "Verify pathway"], pathwayRows))}
        ${detailDisclosure("Associated CpG sites", tableHtml(["CpG site", "Verify"], cpgRows))}
      </section>
    `,
    "gene",
    activeModule,
  );
  if (focus) {
    requestAnimationFrame(() => {
      document.getElementById(focus)?.scrollIntoView({ block: "start" });
    });
  }
}

function renderCpg(record) {
  const sourceDescriptions = (record.sourceDescriptions ?? [])
    .map(
      (item) => `
        <div class="source-description-item">
          <strong>${escapeHtml(item.source)}</strong>
          <p>${escapeHtml(item.text)}</p>
          <div class="external-links">
            ${sourceLink(item.sourceUrl, `${item.sourceSheet} source table`)}
            ${sourceLink(item.publicationUrl, "Publication")}
          </div>
        </div>`,
    )
    .join("");
  const mappingRows = (record.geneMappings ?? []).map((item) => [
    entityLink("gene", item.symbol, item.symbol, "table-link", "epigenomics"),
    escapeHtml(item.mapping),
    sourceLink(item.sourceUrl, "Source"),
  ]);
  const evidenceRows = (record.agingEvidence ?? []).map((item) => [
    escapeHtml(item.source),
    escapeHtml(item.endpoint),
    escapeHtml(item.model),
    escapeHtml(item.effectType),
    sourcedValue(scientific(item.effect), item.sourceUrl),
    sourcedValue(scientific(item.pValue), item.sourceUrl),
    sourceLink(item.sourceUrl, `${item.sourceSheet} source table`),
    sourceLink(item.publicationUrl, "Publication"),
  ]);
  const locationRows = [
    ["Genome build", formatValue(record.genomeBuild)],
    ["Chromosome", formatValue(record.chromosome)],
    ["Position", formatValue(record.position)],
  ];
  shell(
    `
      ${identityHeader("Epigenomics: CpG site", record.id, record.name, record.url, "Source table")}
      <section class="gene-function">
        <div class="eyebrow">Source description</div>
        ${sourceDescriptions}
      </section>
      <section class="detail-stack">
        ${detailDisclosure("Genomic position", tableHtml(["Field", "Value"], locationRows), true)}
        ${detailDisclosure("Source gene annotations", tableHtml(["Gene", "Mapping", "Verify"], mappingRows), true)}
        ${detailDisclosure("Age and mortality associations", tableHtml(["Source", "Endpoint", "Model", "Effect type", "Effect", "P value", "Source table", "Publication"], evidenceRows), true)}
      </section>
    `,
    "cpg",
    "epigenomics",
  );
}

function renderProtein(record) {
  const functionSources = (record.functionPubmedIds ?? []).map((pmid) =>
    externalLink(`https://pubmed.ncbi.nlm.nih.gov/${pmid}/`, `PMID ${pmid}`),
  );
  const organAgeRows = (record.organAgeEvidence ?? []).map((item) => [
    escapeHtml(item.organ),
    escapeHtml(item.targetName),
    sourcedValue(scientific(item.meanNonzeroCoefficient), item.sourceUrl),
    escapeHtml(item.coefficientDirection),
    sourcedValue(formatValue(item.selectedModels), item.sourceUrl),
    sourceLink(item.sourceUrl, "Published models"),
    sourceLink(item.publicationUrl, "Publication"),
  ]);
  const pathwayRows = (record.pathways ?? []).map((item) => [
    entityLink("pathway", item.id, item.id, "table-link"),
    escapeHtml(item.name),
    sourceLink(item.sourceUrl, "Reactome entity"),
    sourceLink(item.pathwayUrl, "Pathway"),
  ]);
  const geneRows = record.geneSymbol
    ? [[
        entityLink("gene", record.geneSymbol, record.geneSymbol, "table-link"),
        sourceLink(record.url, "UniProtKB"),
      ]]
    : [];
  shell(
    `
      ${identityHeader("Proteomics", record.id, record.name, record.url, "UniProt accession")}
      <section class="gene-function">
        <div class="eyebrow">Protein function</div>
        ${sourceTextBlock(record.description, record.url, record.descriptionSource)}
        <div class="external-links">
          ${functionSources.join(" ")}
        </div>
      </section>
      <section class="detail-stack">
        ${detailDisclosure("Encoding gene", tableHtml(["Gene", "Verify"], geneRows), true)}
        ${detailDisclosure("OrganAge evidence", tableHtml(["Organ", "Target", "Mean nonzero coefficient", "Direction", "Selected models", "Published models", "Publication"], organAgeRows), true)}
        ${detailDisclosure("Direct Reactome pathways", tableHtml(["Reactome ID", "Pathway", "Verify entity", "Verify pathway"], pathwayRows))}
      </section>
    `,
    "protein",
    "proteomics",
  );
}

function renderMetabolite(record) {
  const identityRows = [
    ...(record.formula ? [["Formula", sourcedValue(escapeHtml(record.formula), record.url)]] : []),
    ...(record.averageMass ? [["Average mass", sourcedValue(`${escapeHtml(record.averageMass)} Da`, record.url)]] : []),
    ...(record.monoisotopicMass ? [["Monoisotopic mass", sourcedValue(`${escapeHtml(record.monoisotopicMass)} Da`, record.url)]] : []),
    ...(record.charge ? [["Charge", sourcedValue(escapeHtml(record.charge), record.url)]] : []),
  ];
  const evidenceRows = (record.agingEvidence ?? []).map((item) => [
    formatValue(item.method),
    sourcedValue(formatValue(item.value), item.sourceUrl),
    formatValue(item.uncertainty),
    formatValue(item.unit),
    formatValue(item.ageGroup),
    formatValue(item.sample),
    sourceLink(item.sourceUrl, `${item.sourceSheet} source table`),
    item.pmid
      ? externalLink(item.publicationUrl, `PMID ${item.pmid}`)
      : "—",
  ]);
  const pathwayRows = (record.pathways ?? []).map((item) => [
    entityLink("pathway", item.id, item.id, "table-link"),
    escapeHtml(item.name),
    sourceLink(item.sourceUrl, "Reactome entity"),
    sourceLink(item.pathwayUrl, "Pathway"),
  ]);
  shell(
    `
      ${identityHeader("Metabolomics", record.id, record.name, record.url, "ChEBI ID")}
      <section class="gene-function">
        <div class="eyebrow">Source description</div>
        ${sourceTextBlock(record.description, record.url, record.descriptionSource)}
      </section>
      <section class="detail-stack">
        ${detailDisclosure("ChEBI properties", tableHtml(["Property", "Source value"], identityRows), true)}
        ${detailDisclosure("MetaboAge evidence", tableHtml(["Method", "Reported value", "Uncertainty", "Unit", "Age group", "Sample", "Source table", "Publication"], evidenceRows), true)}
        ${detailDisclosure("Direct Reactome pathways", tableHtml(["Reactome ID", "Pathway", "Verify entity", "Verify pathway"], pathwayRows))}
      </section>
    `,
    "metabolite",
    "metabolomics",
  );
}

function renderPathway(record) {
  const geneRows = (record.genes ?? []).map((item) => [
    entityLink("gene", item.symbol, item.symbol, "table-link"),
    escapeHtml(item.name),
    escapeHtml(item.physicalEntity),
    evidenceBadges(item.evidenceLayers, true),
    sourceLink(item.sourceUrl, "Reactome entity"),
  ]);
  const proteinRows = (record.proteins ?? []).map((item) => [
    entityLink("protein", item.id, item.id, "table-link"),
    escapeHtml(item.name),
    sourceLink(item.sourceUrl, "Reactome entity"),
  ]);
  const metaboliteRows = (record.metabolites ?? []).map((item) => [
    entityLink("metabolite", item.id, item.id, "table-link"),
    escapeHtml(item.name),
    sourceLink(item.sourceUrl, "Reactome entity"),
  ]);
  shell(
    `
      ${identityHeader("Pathways", record.id, record.name, record.recordUrl, "Reactome ID")}
      <section class="gene-function">
        <div class="eyebrow">Reactome description</div>
        ${sourceTextBlock(record.description, record.recordUrl, record.descriptionSource)}
      </section>
      <section class="detail-stack">
        ${detailDisclosure("Directly mapped genes", tableHtml(["Gene", "Approved name", "Reactome entity", "Aging evidence", "Verify"], geneRows), true)}
        ${detailDisclosure("Directly mapped proteins", tableHtml(["UniProt accession", "Reactome entity", "Verify"], proteinRows))}
        ${detailDisclosure("Directly mapped metabolites", tableHtml(["ChEBI ID", "Reactome entity", "Verify"], metaboliteRows))}
      </section>
    `,
    "pathway",
    "pathways",
  );
}

function connectionsMap() {
  return `
    <section class="connection-panel atlas-connection">
      <div class="section-heading connection-heading">
        <div><div class="eyebrow">Atlas structure</div><h2>Connected evidence</h2></div>
        <p>Select an area to inspect its source records.</p>
      </div>
      <div class="connection-map" aria-label="Connected human evidence map">
        <svg viewBox="0 0 1260 390" role="img" aria-label="Gene connects to transcript and protein. Epigenomics annotates the transition from gene to transcript. Protein and metabolite records connect directly to pathways.">
          <defs>
            <marker id="arrow" markerWidth="5" markerHeight="5" refX="4.5" refY="2.5" orient="auto">
              <path d="M0,0 L5,2.5 L0,5 Z"></path>
            </marker>
          </defs>
          <rect class="map-field" x="20" y="20" width="1220" height="350" rx="5"></rect>
          <rect class="dogma-field" x="48" y="55" width="770" height="190" rx="5"></rect>
          <text class="map-section-label" x="72" y="82">CENTRAL DOGMA</text>
          <line class="map-section-rule" x1="72" y1="96" x2="794" y2="96"></line>

          <line class="map-arrow" x1="250" y1="165" x2="365" y2="165"></line>
          <line class="map-arrow" x1="525" y1="165" x2="640" y2="165"></line>
          <path class="map-arrow" d="M250 310 C300 310 300 205 337 166"></path>
          <line class="map-arrow" x1="800" y1="165" x2="930" y2="165"></line>
          <path class="map-arrow" d="M800 310 C865 310 875 205 930 172"></path>

          <a href="#/module/genomics" class="svg-node" aria-label="Open Genomics">
            <rect x="90" y="126" width="160" height="78"></rect>
            <line class="node-accent gold" x1="90" y1="126" x2="250" y2="126"></line>
            <text class="node-title" x="170" y="153">Gene</text>
            <text class="node-source node-source-compact" x="170" y="173">
              <tspan x="170" dy="0">GenAge</tspan>
              <tspan x="170" dy="14">LongevityMap</tspan>
            </text>
          </a>
          <a href="#/module/transcriptomics" class="svg-node" aria-label="Open Transcriptomics">
            <rect x="365" y="126" width="160" height="78"></rect>
            <line class="node-accent teal" x1="365" y1="126" x2="525" y2="126"></line>
            <text class="node-title" x="445" y="153">Transcript</text>
            <text class="node-source" x="445" y="177">Human tAge</text>
          </a>
          <a href="#/module/proteomics" class="svg-node" aria-label="Open Proteomics">
            <rect x="640" y="126" width="160" height="78"></rect>
            <line class="node-accent teal" x1="640" y1="126" x2="800" y2="126"></line>
            <text class="node-title" x="720" y="153">Protein</text>
            <text class="node-source" x="720" y="177">OrganAge</text>
          </a>
          <a href="#/module/epigenomics" class="svg-node" aria-label="Open Epigenomics">
            <rect x="90" y="271" width="160" height="78"></rect>
            <line class="node-accent teal" x1="90" y1="271" x2="250" y2="271"></line>
            <text class="node-title" x="170" y="298">Epigenomics</text>
            <text class="node-source" x="170" y="322">cAge · bAge</text>
          </a>
          <a href="#/module/metabolomics" class="svg-node" aria-label="Open Metabolomics">
            <rect x="640" y="271" width="160" height="78"></rect>
            <line class="node-accent teal" x1="640" y1="271" x2="800" y2="271"></line>
            <text class="node-title" x="720" y="298">Metabolite</text>
            <text class="node-source" x="720" y="322">MetaboAge DB</text>
          </a>
          <a href="#/module/pathways" class="svg-node" aria-label="Open Pathways">
            <rect x="930" y="126" width="160" height="78"></rect>
            <line class="node-accent teal" x1="930" y1="126" x2="1090" y2="126"></line>
            <text class="node-title" x="1010" y="153">Pathway</text>
            <text class="node-source" x="1010" y="177">Reactome</text>
          </a>
        </svg>
      </div>
      <p class="connection-note">Arrows show the direct source mappings used for navigation. They do not claim biological effect or causality.</p>
    </section>`;
}

function renderMethods() {
  const sourceByName = Object.fromEntries(
    state.collections.sources.map((source) => [source.name, source]),
  );
  const rows = state.collections.modules.map((module) => [
    escapeHtml(module.label),
    module.sources
      .map((name) => externalLink(sourceByName[name]?.url, name))
      .join(" + "),
    escapeHtml(module.description),
  ]);
  shell(
    `
      ${pageHeader(
        "Scope and interpretation",
        "Methods",
        "A reduced static reference in which every displayed scientific statement is tied to a retained source.",
      )}
      <section class="prose-section">
        <h2>Gene selection</h2>
        <p>A fixed human gene panel retains established aging and longevity searches. Genes are selected reproducibly by evidence-layer breadth, source breadth, source-record depth, and gene symbol. No evidence values are entered by hand.</p>
      </section>
      <section class="prose-section">
        <h2>Published areas</h2>
        ${tableHtml(["Area", "Retained source", "Use"], rows, { limit: 20 })}
      </section>
      <section class="prose-section">
        <h2>Direct pathway mappings</h2>
        <p>${externalLink("https://reactome.org/", "Reactome")} supplies the human gene, protein, and metabolite pathway mappings. The atlas retains direct physical-entity mappings only. Shared pathway membership is not converted into a gene-metabolite relationship.</p>
      </section>
      <section class="prose-section">
        <h2>Scientific boundaries</h2>
        <ul class="method-list">
          ${state.manifest.scientificRules.map((rule) => `<li>${escapeHtml(rule)}</li>`).join("")}
        </ul>
      </section>
    `,
    "methods",
  );
}

function renderSources() {
  const groups = state.collections.sources.reduce((result, item) => {
    (result[item.group] ??= []).push(item);
    return result;
  }, {});
  shell(
    `
      ${pageHeader(
        "Source landscape",
        "Sources",
        "Evidence sources and the identity resources required to interpret them.",
      )}
      ${Object.entries(groups)
        .map(
          ([group, sources]) => `
            <section class="source-group">
              <div class="eyebrow">${escapeHtml(group)}</div>
              ${sources
                .map(
                  (source) => `
                    <div class="source-row">
                      <div><h2>${escapeHtml(source.name)}</h2><p>${escapeHtml(source.role)}</p></div>
                      <div class="source-actions">
                        ${externalLink(
                          source.url,
                          source.name === "tAge"
                            ? "Supplementary Table 2 (XLSX)"
                            : "Open data source",
                        )}
                        ${source.publicationUrl ? externalLink(source.publicationUrl, "Publication") : ""}
                      </div>
                    </div>`,
                )
                .join("")}
            </section>`,
        )
        .join("")}
    `,
    "sources",
  );
}

function renderNotFound(type, id) {
  shell(
    `
      ${pageHeader(
        "Record not found",
        id,
        `No indexed ${ENTITY_LABELS[type]?.toLowerCase() ?? "record"} was found under this identifier.`,
      )}
      <a class="primary-link" href="#/atlas">Return to Atlas search</a>
    `,
    type,
    ENTITY_MODULES[type] ?? "",
  );
}

async function renderEntity(type, id, focus = "") {
  const activeModule =
    focus && moduleDefinition(focus) ? focus : ENTITY_MODULES[type];
  shell(
    `<div class="loading-state">Loading ${escapeHtml(ENTITY_LABELS[type])}...</div>`,
    type,
    activeModule,
  );
  const record = await loadEntity(type, id);
  if (!record) {
    renderNotFound(type, id);
    return;
  }
  if (type === "gene") renderGene(record, focus);
  else if (type === "cpg") renderCpg(record);
  else if (type === "protein") renderProtein(record);
  else if (type === "metabolite") renderMetabolite(record);
  else if (type === "pathway") renderPathway(record);
}

async function route() {
  const { section, id, focus } = parseRoute();
  window.scrollTo(0, 0);
  if (section === "atlas" || section === "connections") renderAtlas();
  else if (section === "module") renderModule(id);
  else if (section === "methods") renderMethods();
  else if (section === "sources") renderSources();
  else if (ENTITY_LABELS[section] && id) await renderEntity(section, id, focus);
  else renderNotFound("record", id || section);
}

async function initialize() {
  try {
    [state.manifest, state.collections, state.searchIndex] = await Promise.all([
      getJson("data/manifest.json"),
      getJson("data/collections.json"),
      getJson("data/search-index.json"),
    ]);
    await route();
    window.addEventListener("hashchange", route);
  } catch (error) {
    app.innerHTML = `
      <div class="fatal-error">
        <h1>Unable to load the atlas</h1>
        <p>${escapeHtml(error.message)}</p>
      </div>`;
  }
}

initialize();
