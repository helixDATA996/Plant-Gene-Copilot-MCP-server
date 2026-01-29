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
Create a virtual environment with `uv` and install `requirements.txt`

```bash
# 1) Install uv
# macOS / Linux
curl -LsSf https://astral.sh/uv/install.sh | sh
# Windows (PowerShell)
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
#Using pip
pip install uv

# 2) Go to your project directory
cd /path/to/your/project

# 3) Create a virtual environment (defaults to .venv)
uv venv

# 4) Activate the environment
# macOS / Linux
source .venv/bin/activate
# Windows (PowerShell)
.venv\Scripts\Activate.ps1

# 5) Install dependencies
uv pip install -r requirements.txt

# 6) Verify
python -V
python -c "import sys; print(sys.executable)"

# 7)Configration
Config your NCBI api key

# 8)Run
python mcp_server.py
