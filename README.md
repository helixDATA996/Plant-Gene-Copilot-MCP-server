# Plant-Gene-Copilot-MCP-server

A lightweight, open-source MCP (Model Context Protocol) frontend for structured plant genomics research.  
It connects to a backend MCP server that provides unified function-style tools for gene lookup, identifier normalization, multi-source annotation retrieval (NCBI/UniProt/KEGG/InterPro/Pfam), genome localization (GFF-derived), and literature search (PubMed), enabling LLM agents to access plant genomic knowledge in a consistent and reproducible way.

## Key Features

- Unified MCP tool interfaces for plant functional genomics workflows
- Gene-centric retrieval: gene search, ID normalization, cross-database aliases
- Annotation integration: NCBI Gene, UniProt, KEGG, InterPro + Pfam/CDD/SMART mappings
- Genome context: chromosome coordinates, neighborhood genes, genomic location visualization
- Protein interaction utilities via STRING (resolution, network, enrichment)
- Literature retrieval via PubMed (basic + advanced query + article details)
- Designed for tool-augmented LLM agents and reproducible bioinformatics pipelines

## Architecture

This repository provides the MCP frontend (MCP tool client).  
The backend MCP server (not included unless you add it) bridges MCP tool calls to an intermediary API service and underlying MySQL repositories.

### High-level flow

1. MCP Frontend (this repo)
2. MCP Server endpoints (tool definitions)
3. Intermediary API service (query optimization + data integrity)
4. MySQL knowledge bases + external sources integrations

### Backend Knowledge Bases (for reference)

- **Functional Knowledgebase (9 tables; 1,535,469 records)**
  - `genes`, `main_entry`, `uniprot`, `kegg_details`, `ncbi_gene_details`, `interpro_main`, `interpro_pfam`, `interpro_cdd`, `interpro_smart`
- **Annotation Module / GFF database (7 species tables; 294,808 records)**
  - GFF-derived gene features for: *Arabidopsis thaliana*, *Zea mays*, *Solanum lycopersicum*, *Populus nigra*, *Populus trichocarpa*, *Oryza sativa ssp. japonica*, *Marchantia polymorpha*
- **Total:** 16 tables; 1,830,277 records (~1.83M)

## Tool List (MCP Functions)

The tools below are grouped by function. Names shown in `code` are the
friendly/logical names; the exact names registered with the MCP client are the
Python function names in `mcp_server.py`, which use `snake_case`
(e.g. `geneSearch` below is exposed as `gene_search`). Verify against the
`@mcp.tool()` decorated functions in the source before scripting.

### 1) Gene lookup & identifier normalization

- **`geneSearch(gene_name, extended_search=false)`**  
  Search genes by symbol/name/gene_id; supports multiple queries separated by `;`.  
  Returns curated gene records and RefSeq links when available.

- **`fetchGeneId(gene_name)`**  
  Resolve gene name → gene_id(s).

- **`entrezgeneId(gene_id)`**  
  Convert internal gene_id (e.g., `AT1G01060`) → NCBI Entrez Gene numeric ID (required by NCBI tools).

- **`fetchExternaldbs(gene_id, max_results=10, compact=true)`**  
  Retrieve aliases/synonyms and cross-database identifiers (TAIR/Araport/etc.) from NCBI mappings.

### 2) Genome annotation & localization (GFF)

- **`fetchGffAnnotation(gene_id)`**  
  Returns chromosome coordinates, genomic neighborhood context, and genome-location visualization URL(s).

### 3) NCBI Gene: structure & transcripts

- **`fetchNcbiInfo(entrezgene_id)`**  
  Gene summary and an exon–intron structure visualization (image link) from NCBI Gene.

- **`fetchGeneStructureDetails(entrezgene_id)`**  
  Transcript/exon/intron coordinate details for alternative splicing analyses.

### 4) UniProt protein metadata

- **`fetchUniprotEntry(gene_id, extentional_search=false)`**  
  Resolve gene_id(s) → UniProt entry/entries.

- **`fetchUniprotDetails(entries)`**  
  Retrieve detailed UniProt annotations for entries (includes UniProt links).

### 5) KEGG functional annotation

- **`fetchKeggInfo(kegg_id)`**  
  Retrieve pathways, orthology, GO terms, enzyme/reaction, and DB links.

### 6) InterPro / Pfam domain & family annotation

- **`fetchInterpro(interpro)`**  
  InterPro entry details.

- **`fetchPfam(pfam)`**  
  Pfam domain details.

### 7) PubMed literature retrieval

- **`fetchPmid(query, max_results=20)`**  
  Basic PubMed query → PMIDs and titles.

- **`pubmedAdvancedSearch(query, page=1, sort="relevance", date_from="", date_to="", article_type="")`**  
  Advanced PubMed search with pagination and filters.

- **`fetchDetails(pmid)`**  
  Retrieve full article metadata (title/authors/abstract/DOI/etc.) for PMID list.

### 8) STRING protein interaction utilities

- **`stringQuerySpecies(species_text)`**  
  Resolve species name/taxon hints for STRING.

- **`stringResolveProteins(proteins, species, show_sequence)`**  
  Map identifiers → STRING IDs.

- **`stringInteractionsQuerySet(proteins, species, required_score, network_type, extend_network)`**  
  Get interactions among a given query set.

- **`stringAllInteractionPartners(identifiers, species, required_score, network_type)`**  
  Retrieve all partners for proteins.

- **`stringEnrichment(proteins, species)`**  
  Functional enrichment over a protein list.

- **`stringFunctionalAnnotation(identifiers, species)`**  
  Retrieve functional annotations.

- **`stringVisualNetwork(proteins, species, required_score, network_type, network_flavor, extend_network, hide_disconnected_nodes)`**  
  Generate a network image URL.

- **`stringNetworkLink(proteins, species, required_score, network_type, network_flavor, extend_network, hide_disconnected_nodes)`**  
  Generate an interactive network link.

- **`stringInteractionEvidence(identifier_a, identifiers_b, species)`**  
  Evidence pages for edges.

- **`stringHomology(proteins, species, species_b)`**  
  Homolog lookup.

- **`stringPpiEnrichment(identifiers, species, required_score)`**  
  Check PPI enrichment.

- **`stringProteinsForTerm(species, term_text)`**  
  Proteins annotated with a given functional term.

- **`stringSequenceSearch(sequences, species)`**  
  Sequence-based lookup in STRING.

- **`stringEnrichmentImageUrl(identifiers, species, category)`**  
  Enrichment plot URL.

- **`stringHelp(topic)`**  
  Help for STRING tool features.

### 9) Service health

- **`checkInitialization()`**  
  Whether the backend services are initialized and reachable.

## Quick Start

The server needs **Python 3.10+**. It talks to an existing backend API
(`http://api.plantgenecopilot.top`, overridable in `config/server.config`) plus
the public STRING and NCBI E-utilities services.

### 1. Install dependencies

Using `uv` (recommended):

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh        # macOS / Linux
# or: powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

Or with `pip` + `venv` (no extra tooling):

```bash
python -m venv .venv
source .venv/bin/activate          # macOS / Linux
# .venv\Scripts\Activate.ps1        # Windows (PowerShell)
pip install -r requirements.txt
```

### 2. Configure API access

Edit `config/server.config` before running:

- `ncbi_config.email` — a valid email (required by NCBI E-utilities policy).
- `ncbi_config.api_key` — your NCBI API key (optional, raises rate limits).
  Obtainable at <https://www.ncbi.nlm.nih.gov/account/settings/>.
- `api_base_url` — backend API address. Change it if you run your own backend
  instead of the hosted one. The `base_url`, `server_port` and `timeout`
  entries tune the STRING service.

### 3. Run

```bash
python mcp_server.py
```

The server starts over **HTTP transport** on `0.0.0.0:1154` (see the log
output). Point any MCP-compatible client (Claude Desktop, Dify, or a FastMCP
client) at that endpoint to begin querying plant-genomics tools. Press `q`
to stop.

### Calling a tool from your own code

With a FastMCP/anyio MCP client, connect to the endpoint and call a tool by
its registered name, for example `gene_search(gene_name="AtDCL1")`:

```python
import asyncio
from fastmcp import Client

async def main():
    async with Client("http://localhost:1154") as client:
        resp = await client.call_tool("gene_search", {"gene_name": "AtDCL1"})
        print(resp)

asyncio.run(main())
```

> Note: tool names use `snake_case` as registered in `mcp_server.py`
> (e.g. `gene_search`, `fetch_gene_id`, `fetch_pmid`). The full annotated
> tool list is in the next section.


<p>
  本项目的 AI API 支持由
  <a href="https://tokeness.io">
    Tokeness.io
  </a>
  赞助提供。
</p>


## Reproducibility

The two evaluation components of the manuscript are separated into
`workflow_evaluation/` (automated multi-choice benchmark) and
`expert_evaluation/` (human expert panel).

### Automated workflow evaluation — `workflow_evaluation/`

Machine-readable benchmark and per-question outputs, fully reproducible:

- `datasets/` — the Expert MoBiPlant benchmark as JSON:
  `expert_mobi.json` (full 565 questions), `filtered_species.json`
  (308-question target-plant subset), `other_species.json` (257-question
  non-target subset). Each entry: `question`, `options` (3 choices),
  `answer`, `plant_species`, `area` (domain), `source`/`doi`.
- `model_outputs/` — per-question outputs for 6 models x 2 datasets
  (agent mode = LLM + MCP tools; LLM mode = no tools), with judge and
  objective correctness.
- `summary_tables/` — Table 1 (performance comparison with McNemar
  significance), Table 2 (area breakdown), species-stratified accuracies
  with 95% Wilson CI, per-domain sample sizes, and a per-question answer
  table.
- `analysis/` — scripts: `run_dify_evaluation.py` (runs the benchmark via
  the Dify workflow API), `eval_mcp_workflow_target.yml` (exported Dify
  workflow DSL, contains no credentials), `compute_summary_tables.py`
  (recomputes the summary tables from `model_outputs/`), and `make_fig2.py`
  / `make_fig3.py`.

### Human expert evaluation — `expert_evaluation/`

- `system/` — the web-based double-blind scoring system used by the 15
  expert assessors (Flask + MySQL).
- `data/expert_scores_anonymized.xlsx` — the full set of 1,050 anonymized
  ratings (15 raters x 10 questions x 7 metrics), with the ICC,
  significance-test, and rubric sheets.
- `analysis/mixed_effects_bootstrap.py` — mixed-effects and clustered
  bootstrap analyses of the expert ratings.

A versioned snapshot of this repository is archived on Zenodo (DOI in the
manuscript).
