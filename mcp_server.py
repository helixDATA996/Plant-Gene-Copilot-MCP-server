"""
Plant Gene Database MCP Server
==============================

Lightweight MCP server that provides tools for querying plant gene information.
All database operations are delegated to the HTTP API server (http_server.py).

Run:
  python mcp_server.py

Requires:
  - http_server.py running on port 8080
  - fastmcp>=2.10.6
  - httpx>=0.28
"""

import sys
import json
import time
import httpx
import asyncio
import logging
import threading
import os
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Annotated, Optional
from collections import defaultdict
from pydantic import Field
from fastmcp import FastMCP
from fastmcp.server.dependencies import get_http_headers

# ============================================================================
# Logging Configuration
# ============================================================================

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stderr)]
)
logger = logging.getLogger(__name__)

# ============================================================================
# Configuration
# ============================================================================

CONFIG_FILE = "config/server.config"

# HTTP API Server URL (can be overridden by config file or env var)
API_BASE_URL = "http://localhost:8080"

# STRING API Configuration (loaded from config)
string_base_url = None
string_timeout = 100.0
log_verbosity = {'call': False, 'params': False, 'taskid': False, 'size': False}

# NCBI API Configuration (loaded from config, for PubMed search)p
ncbi_email = None
ncbi_api_key = None
ncbi_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def load_configuration():
    """Load configuration from config file"""
    global string_base_url, string_timeout, log_verbosity, API_BASE_URL
    global ncbi_email, ncbi_api_key, ncbi_base_url

    config_path = Path(CONFIG_FILE)
    if not config_path.exists():
        logger.warning(f"Config file not found: {CONFIG_FILE}")
        return False

    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            config = json.load(f)

        # Load API base URL (env var takes priority)
        if os.getenv("API_BASE_URL"):
            API_BASE_URL = os.getenv("API_BASE_URL")
        elif 'api_base_url' in config:
            API_BASE_URL = config.get("api_base_url")
        logger.info(f"API Base URL: {API_BASE_URL}")

        # Load STRING API configuration
        if 'base_url' in config:
            string_base_url = config.get("base_url")
            string_timeout = float(config.get("timeout", 100))
            if 'verbosity' in config:
                if config['verbosity'] == 'full':
                    log_verbosity = {'call': True, 'params': True, 'taskid': True, 'size': True}
                elif config['verbosity'] == 'low':
                    log_verbosity = {'call': True, 'params': False, 'taskid': True, 'size': True}
            logger.info(f"STRING API configured: {string_base_url}")

        # Load NCBI API configuration (env vars take priority)
        if 'ncbi_config' in config:
            ncbi_cfg = config['ncbi_config']
            ncbi_email = os.getenv("NCBI_EMAIL") or ncbi_cfg.get('email')
            ncbi_api_key = os.getenv("NCBI_API_KEY") or ncbi_cfg.get('api_key')
            ncbi_base_url = ncbi_cfg.get('base_url', ncbi_base_url)
            if ncbi_email and ncbi_api_key:
                logger.info(f"NCBI API configured: {ncbi_email}")
            else:
                logger.warning("NCBI API not fully configured (missing email or api_key)")
        else:
            # Try env vars only
            ncbi_email = os.getenv("NCBI_EMAIL")
            ncbi_api_key = os.getenv("NCBI_API_KEY")
            if ncbi_email and ncbi_api_key:
                logger.info(f"NCBI API configured from env: {ncbi_email}")

        return True
    except Exception as e:
        logger.error(f"Error loading config: {e}")
        return False


load_configuration()

# ============================================================================
# STRING Help Topics
# ============================================================================

HELP_TOPICS = {
    "gsea": "GSEA cannot be performed directly by the agent. Use the 'Proteins with Values/Ranks' option on the STRING input page with the complete set of proteins from your experiment.",
    "clustering": "Clustering is available on the STRING website under the Clustering tab. STRING offers MCL and k-means methods.",
    "large_input": "Large input sets may cause timeouts. Use the STRING web interface or Cytoscape STRING app for large-scale analyses.",
    "scores": "STRING interaction scores range from 0-1000. Common thresholds: 400=medium, 700=high confidence.",
    "missing_proteins": "STRING accepts many identifiers. If a protein cannot be found, check the older STRING v11.5.",
    "missing_species": "Use the 'Add species' functionality on STRING to upload custom proteomes.",
    "proteome_annotation": "Use the 'Add species' functionality for proteome annotation with GO and KEGG predictions.",
    "regulatory_networks": "Regulatory/directed networks are not available in STRING. All links are undirected.",
    "how_to_use_string": "STRING explores protein-protein interactions and functional enrichment. Provide proteins to retrieve interactions and enrichment results.",
    "line_colors": "STRING networks use confidence view (line thickness) or evidence view (colored by source type).",
}

# ============================================================================
# HTTP API Helper - Global Async Client with Connection Pool
# ============================================================================

# Global async HTTP client for internal API (connection pool reuse)
_http_client: httpx.AsyncClient = None
_http_client_lock: asyncio.Lock = None

# Global async HTTP client for STRING API (connection pool reuse)
_string_client: httpx.AsyncClient = None
_string_client_lock: asyncio.Lock = None


def _ensure_locks_initialized():
    """Ensure locks are initialized (call from event loop context)"""
    global _http_client_lock, _string_client_lock
    if _http_client_lock is None:
        _http_client_lock = asyncio.Lock()
    if _string_client_lock is None:
        _string_client_lock = asyncio.Lock()


def _get_http_client_lock() -> asyncio.Lock:
    """Get lock for HTTP client"""
    global _http_client_lock
    if _http_client_lock is None:
        _http_client_lock = asyncio.Lock()
    return _http_client_lock


def _get_string_client_lock() -> asyncio.Lock:
    """Get lock for STRING client"""
    global _string_client_lock
    if _string_client_lock is None:
        _string_client_lock = asyncio.Lock()
    return _string_client_lock


async def get_http_client() -> httpx.AsyncClient:
    """Get or create global HTTP client with connection pool for internal API"""
    global _http_client

    # Fast path: client exists and is open
    if _http_client is not None and not _http_client.is_closed:
        return _http_client

    # Slow path: need to create client (with lock to prevent race condition)
    async with _get_http_client_lock():
        # Double-check after acquiring lock
        if _http_client is not None and not _http_client.is_closed:
            return _http_client

        # Close old client if exists but broken
        if _http_client is not None:
            try:
                await _http_client.aclose()
            except Exception:
                pass

        _http_client = httpx.AsyncClient(
            base_url=API_BASE_URL,
            timeout=httpx.Timeout(120.0, connect=15.0),  # 总超时120秒，连接超时15秒
            limits=httpx.Limits(
                max_connections=200,
                max_keepalive_connections=50,
                keepalive_expiry=30.0
            )
        )
        logger.info(f"HTTP client created: {API_BASE_URL}")
    return _http_client


async def get_string_client() -> httpx.AsyncClient:
    """Get or create global HTTP client with connection pool for STRING API"""
    global _string_client

    if not string_base_url:
        raise RuntimeError("STRING API not configured")

    # Fast path
    if _string_client is not None and not _string_client.is_closed:
        return _string_client

    # Slow path with lock
    async with _get_string_client_lock():
        if _string_client is not None and not _string_client.is_closed:
            return _string_client

        if _string_client is not None:
            try:
                await _string_client.aclose()
            except Exception:
                pass

        _string_client = httpx.AsyncClient(
            base_url=string_base_url,
            timeout=httpx.Timeout(string_timeout, connect=15.0),  # 连接超时15秒
            limits=httpx.Limits(
                max_connections=100,
                max_keepalive_connections=20,
                keepalive_expiry=60.0
            )
        )
        logger.info(f"STRING client created: {string_base_url}")
    return _string_client


async def call_api(endpoint: str, params: dict = None, retries: int = 2) -> dict:
    """Async API call with retry logic"""
    last_error = None

    for attempt in range(retries + 1):
        try:
            client = await get_http_client()
            resp = await client.get(f"/api/{endpoint}", params=params)
            resp.raise_for_status()
            return resp.json()
        except httpx.TimeoutException as e:
            logger.error(f"API call timeout: {endpoint} - {e}")
            return {"error": f"请求超时: {endpoint} (超过120秒)，请稍后重试或减少查询数量"}
        except (httpx.ConnectError, httpx.PoolTimeout, httpx.RemoteProtocolError) as e:
            last_error = e
            logger.warning(f"API call attempt {attempt + 1} failed: {endpoint} - {e}")
            # Reset client on connection errors
            await reset_http_client()
            if attempt < retries:
                await asyncio.sleep(0.5 * (attempt + 1))  # Exponential backoff
        except Exception as e:
            logger.error(f"API call failed: {endpoint} - {e}")
            return {"error": str(e)}

    logger.error(f"API call failed after {retries + 1} attempts: {endpoint} - {last_error}")
    return {"error": f"连接失败: {last_error}"}


async def reset_http_client():
    """Reset HTTP client (force reconnection)"""
    global _http_client
    async with _get_http_client_lock():
        if _http_client is not None:
            try:
                await _http_client.aclose()
            except Exception:
                pass
            _http_client = None
            logger.info("HTTP client reset")


async def reset_string_client():
    """Reset STRING client (force reconnection)"""
    global _string_client
    async with _get_string_client_lock():
        if _string_client is not None:
            try:
                await _string_client.aclose()
            except Exception:
                pass
            _string_client = None
            logger.info("STRING client reset")


async def cleanup():
    """Cleanup resources on shutdown"""
    global _http_client, _string_client

    if _http_client is not None:
        try:
            await _http_client.aclose()
            logger.info("HTTP client closed")
        except Exception as e:
            logger.warning(f"Error closing HTTP client: {e}")
        _http_client = None

    if _string_client is not None:
        try:
            await _string_client.aclose()
            logger.info("STRING client closed")
        except Exception as e:
            logger.warning(f"Error closing STRING client: {e}")
        _string_client = None


# ============================================================================
# STRING API Helper Functions
# ============================================================================

async def _post_json(client: httpx.AsyncClient, endpoint: str, data: dict, retries: int = 1):
    """POST to STRING API endpoint with retry logic"""
    last_error = None

    for attempt in range(retries + 1):
        try:
            response = await client.post(endpoint, data=data)
            response.raise_for_status()
            try:
                return response.json()
            except ValueError:
                return {"result": response.text}
        except httpx.TimeoutException as e:
            logger.error(f"STRING API timeout: {endpoint} - {e}")
            return {"error": {"type": "timeout_error", "message": f"STRING API请求超时 ({endpoint})，请稍后重试或减少蛋白质数量"}}
        except (httpx.ConnectError, httpx.PoolTimeout, httpx.RemoteProtocolError) as e:
            last_error = e
            logger.warning(f"STRING API attempt {attempt + 1} failed: {endpoint} - {e}")
            if attempt < retries:
                await reset_string_client()
                await asyncio.sleep(0.5 * (attempt + 1))
        except httpx.HTTPStatusError as e:
            return {"error": {"type": "string_api_error", "message": f"STRING API错误: HTTP {e.response.status_code}"}}
        except Exception as e:
            return {"error": {"type": "unexpected_error", "message": str(e)}}

    return {"error": {"type": "connection_error", "message": f"STRING API连接失败: {last_error}"}}


def truncate_enrichment(data, is_json):
    """Truncate enrichment results"""
    if is_json.lower() != 'json':
        return data
    term_cutoff, size_cutoff = 20, 15
    filtered = []
    category_count = defaultdict(int)
    for row in data:
        category = row.get('category', '')
        category_count[category] += 1
        if category_count[category] > term_cutoff:
            continue
        row['inputGenes_total'] = len(row.get('inputGenes', []))
        row['preferredNames_total'] = len(row.get('preferredNames', []))
        if len(row.get('inputGenes', [])) > size_cutoff:
            row['inputGenes'] = row['inputGenes'][:size_cutoff] + ["..."]
            row['preferredNames'] = row['preferredNames'][:size_cutoff] + ["..."]
            row['truncated'] = True
        else:
            row['truncated'] = False
        filtered.append(row)
    return filtered


def truncate_network(data, score_threshold=None, size_cutoff=100, is_json="json"):
    """Truncate network results"""
    if is_json.lower() != "json":
        return data, None
    try:
        threshold = float(score_threshold) if score_threshold else 400.0
    except:
        threshold = 400.0
    if threshold > 1:
        threshold = threshold / 1000.0
    filtered = [row for row in data if row.get("score", 0) >= threshold]
    filtered.sort(key=lambda r: r.get("score", 0), reverse=True)
    add_note = len(filtered) > size_cutoff
    if add_note:
        filtered = filtered[:size_cutoff]
    return filtered, add_note


def log_call(endpoint, params):
    """Log API call"""
    if log_verbosity['call']:
        logger.info(f"Call: {endpoint}")


def log_response_size(resp):
    """Log response size"""
    pass


def object_size(obj):
    """Calculate object size"""
    if isinstance(obj, str):
        return len(obj)
    elif isinstance(obj, (int, float)):
        return len(str(obj))
    elif isinstance(obj, dict):
        return sum(object_size(v) for v in obj.values())
    elif isinstance(obj, list):
        return sum(object_size(v) for v in obj)
    return 0


# ============================================================================
# Initialize FastMCP Server
# ============================================================================

mcp = FastMCP("Plant Gene & STRING Database MCP Server")

# ============================================================================
# Plant Gene Database Tools (via HTTP API)
# ============================================================================

@mcp.tool()
async def gene_search(gene_name: str, extended_search: bool = False) -> dict:
    """
    Search for plant gene information by gene name (step 1: basic gene lookup)

    This tool automatically tries exact search first, then falls back to fuzzy search if needed.

    Args:
        gene_name:  Gene name, gene symbol or gene_id to search for.Supports semicolon-separated multiple queries (e.g. "ERF4;AT1G01060").
        extended_search: True = allow fuzzy LIKE search if exact search returns empty (default: False)

    Returns:
        Dictionary with:
        - result: List of gene information dictionaries with RefSeq links added (if <=20 results)
          * refseq_mrna_link: NCBI Nucleotide link for mRNA sequence (NM_ accession)
          * refseq_peptide_link: NCBI Protein link for peptide sequence (NP_ accession)
        - result_count: Number of results found
        - search_type: "exact" or "extended" - indicates which search method was used
        - warning: Message if too many results found (>20)
        - note: Instruction for AI model
        - sample_genes: List of 10 sample records with gene_id and gene_names (if >20 results)
    """
    # Step 1: Always try exact search first
    result = await call_api("gene_search", {"gene_name": gene_name})
    search_type = "exact"

    # Handle error responses
    if isinstance(result, dict) and "error" in result:
        return result

    # Step 2: If no results and extended_search is enabled, try fuzzy search
    result_count = len(result) if isinstance(result, list) else 0
    if result_count == 0 and extended_search:
        logger.info(f"No exact match found for '{gene_name}', trying extended search...")
        result = await call_api("extend_api_gene_search", {"gene_name": gene_name})
        search_type = "extended"

        # Handle error responses from extended search
        if isinstance(result, dict) and "error" in result:
            return result

        result_count = len(result) if isinstance(result, list) else 0

    # If too many results (>20), return warning
    if result_count > 20:
        return {
            "warning": f"Too many results found: {result_count} genes match the search criteria.",
            "note": "Too many results. Review the sample_genes below to identify the target gene, then ask the user to provide a more specific gene_id (recommended) or gene name.",
            "result_count": result_count,
            "search_type": search_type,
            "sample_genes": [
                {
                    "gene_id": r.get('Gene_ID', 'N/A'),
                    "gene_names": r.get('Gene_Names', 'N/A')
                }
                for r in result[:10]
            ]
        }

    # If moderate results (<=20), add RefSeq links
    if isinstance(result, list) and result:
        for record in result:
            if isinstance(record, dict):
                # Add RefSeq mRNA link (NM_ accessions)
                refseq_mrna = record.get('refseq_mrna', '')
                if refseq_mrna and refseq_mrna != 'N/A':
                    # Take first RefSeq mRNA ID and create NCBI nucleotide link
                    first_mrna = refseq_mrna.split(',')[0].split(';')[0].strip()
                    if first_mrna:
                        record['refseq_mrna_link'] = f"https://www.ncbi.nlm.nih.gov/nuccore/{first_mrna}"

                # Add RefSeq peptide link (NP_ accessions)
                refseq_peptide = record.get('refseq_peptide', '')
                if refseq_peptide and refseq_peptide != 'N/A':
                    # Take first RefSeq peptide ID and create NCBI protein link
                    first_peptide = refseq_peptide.split(',')[0].split(';')[0].strip()
                    if first_peptide:
                        record['refseq_peptide_link'] = f"https://www.ncbi.nlm.nih.gov/protein/{first_peptide}"

    return {"result": result, "result_count": result_count, "search_type": search_type}


@mcp.tool()
async def fetch_gene_id(gene_name: str) -> str:
    """
    Search by gene name and return semicolon-separated gene_id string

    Args:
        gene_name: Gene name to search, supports semicolon-separated multiple names

    Returns:
        Semicolon-separated gene_id string
    """
    result = await call_api("fetch_gene_id", {"gene_name": gene_name})
    return result.get("result", "")


@mcp.tool()
async def entrezgene_id(gene_id: str) -> str:
    """
    Convert Gene_ID to NCBI Entrezgene ID (required for fetch_ncbi_info)

    Use this tool to get the numeric NCBI ID needed for fetch_ncbi_info.
    Gene_ID (e.g., "AT1G01060") → Entrezgene ID (e.g., "839341")

    Args:
        gene_id: Gene ID from gene_search or database (e.g., "AT1G01060"), supports semicolon-separated multiple IDs

    Returns:
        Semicolon-separated NCBI Entrezgene ID string (e.g., "839341;839342")
        Use these IDs directly with fetch_ncbi_info to get gene structure visualization
    """
    result = await call_api("entrezgene_id", {"gene_id": gene_id})
    return result.get("result", "N/A")

@mcp.tool()
async def fetch_ExternalDBs(
    gene_id: str,
    max_results: int = 10,
    compact: bool = True
) -> dict:
    """
    Find all gene names, aliases and cross-references from NCBI (comprehensive name lookup)

    OPTIMIZED: This tool now limits results and truncates long fields to prevent context overflow.

    Use this tool when you need to find:
    - All alternative names/aliases for a gene
    - Cross-database IDs (TAIR, Araport, etc.)
    - Official nomenclature symbols and names
    - Gene synonyms

    Args:
        gene_id: Gene ID, symbol, or any alias (e.g., "AT1G01060", "DCL1"), supports semicolon-separated IDs
        max_results: Maximum number of results to return (default 10, max 100). Use lower values for broad searches.
        compact: If True (default), truncates long fields (ExternalDBs, Description) to save context tokens.

    Returns:
        Dictionary with:
        - results: List of gene records (limited to max_results)
        - returned: Number of results returned
        - total: Total matching records in database
        - truncated: True if more results exist than returned
        - note: Guidance message if results were truncated
    """
    result = await call_api("fetch_ExternalDBs", {
        "gene_id": gene_id,
        "max_results": max_results,
        "compact": str(compact).lower()
    })
    return result

@mcp.tool()
async def fetch_gff_annotation(gene_id: str) -> dict:
    """
    Get GFF annotation and generate GENOMIC LOCATION visualization (gene position on chromosome)

    Use this for: chromosome position, genomic coordinates, gene neighbors, upstream/downstream genes
    NOT for: exon-intron structure (use fetch_ncbi_info instead)

    Args:
        gene_id: Gene ID (e.g., "AT1G01060"), supports semicolon-separated multiple IDs

    Returns:
        Dictionary with:
        - summary: Formatted markdown summary text
        - result: Clean annotation data (internal fields removed)
        - image_url: List of genomic location diagram URLs
        Picture URLs must be displayed in markdown format if exists
    """
    # Call HTTP API to get raw data
    raw_response = await call_api("fetch_gff_annotation", {"gene_id": gene_id})

    if "error" in raw_response:
        return raw_response

    raw_data = raw_response.get("raw_data", [])
    image_filenames = raw_response.get("image_filenames", [])

    # Separate target and neighbor genes
    target_genes = [r for r in raw_data if r.get('_is_target', True)]

    # Build image URLs (business logic in MCP layer)
    image_urls = [f"{API_BASE_URL}/location_image/{fname}" for fname in image_filenames]

    # Build summary (business logic in MCP layer)
    summary_lines = [f"# GFF Annotation Results for {gene_id}\n",
                     f"**Found {len(target_genes)} annotation record(s)**\n"]

    for idx, rec in enumerate(target_genes, 1):
        start = rec.get('Start', 'N/A')
        end = rec.get('End', 'N/A')
        start_str = f"{start:,}" if isinstance(start, (int, float)) else str(start)
        end_str = f"{end:,}" if isinstance(end, (int, float)) else str(end)
        cds = rec.get('CDS_Length', 0)
        cds_str = f"{cds:,}" if isinstance(cds, (int, float)) else str(cds)

        summary_lines.extend([
            f"\n## Record {idx}",
            f"- **Gene ID**: {rec.get('GeneID', 'N/A')}",
            f"- **Gene Name**: {rec.get('Gene_Name', 'N/A')}",
            f"- **Chromosome**: {rec.get('Chromosome', 'N/A')}",
            f"- **Location**: {start_str} - {end_str} bp ({rec.get('Strand', 'N/A')})",
            f"- **Biotype**: {rec.get('Biotype', 'N/A')}",
            f"- **Source**: {rec.get('Source', 'N/A')}",
            f"- **Transcripts**: {rec.get('Transcript_Count', 0)}",
            f"- **Exons**: {rec.get('Exon_Count', 0)}",
            f"- **Proteins**: {rec.get('Protein_Count', 0)}",
            f"- **CDS Length**: {cds_str} bp",
        ])

        if rec.get('Description'):
            summary_lines.append(f"- **Description**: {rec.get('Description', 'N/A')}")

        # Add table name for reference
        summary_lines.append(f"- **Table**: {rec.get('_table_name', 'N/A')}")

    # Clean result (remove internal fields - business logic)
    clean_targets = [{k: v for k, v in t.items() if not k.startswith('_')} for t in target_genes]

    return {
        "summary": "\n".join(summary_lines),
        "result": clean_targets,
        "image_url": image_urls,
        "notes": ["Gene location visualization generated.", "Embed image links as markdown."]
    }


@mcp.tool()
async def fetch_ncbi_info(entrezgene_id: str) -> list:
    """
    Get NCBI gene info with GENE STRUCTURE visualization (exon-intron architecture)

    Use this for: gene structure, exon-intron organization, transcript variants, CDS regions
    NOT for: genomic position on chromosome (use fetch_gff_annotation instead)

    Args:
        entrezgene_id: NCBI Entrezgene ID (e.g., "839341"), supports semicolon-separated multiple IDs
        Get entrezgene_id from gene_search results or use entrezgene_id() tool to convert Gene_ID

    Returns:
        List of dictionaries with:
        - ncbi_link: NCBI gene database URL
        - structure_image: Gene structure diagram (markdown format, MUST display if exists)
          * Rectangles = Exons (coding regions)
          * V-shaped lines = Introns (spliced out regions)
          * Different colors = Different transcript variants
          * 5' and 3' ends marked with arrows showing transcription direction
        - Gene information: Symbol, Description, Chromosome, Summary, GO terms, etc. (17 key fields)
    """
    # Fields to keep (exclude large fields like ExonDetails, IntronDetails, etc.)
    keep_fields = {
        'entrezgene_id', 'Symbol', 'Description', 'Chromosome', 'NomenclatureSymbol',
        'NomenclatureName', 'OtherAliases', 'Summary', 'MapLocation', 'GeneType',
        'Organism', 'CommonName', 'GenomicStart', 'GenomicStop', 'ExonCount',
        'TranscriptCount', 'GOTerms'
    }

    result = await call_api("fetch_ncbi_info", {"entrezgene_id": entrezgene_id})
    filtered_result = []

    if isinstance(result, list):
        for record in result:
            if isinstance(record, dict):
                # Filter fields
                filtered_record = {k: v for k, v in record.items() if k in keep_fields}
                # Add NCBI link and structure image
                gene_id = record.get('entrezgene_id')
                symbol = record.get('Symbol', '')
                if gene_id:
                    filtered_record['ncbi_link'] = f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}"
                    # Add gene structure image link
                    if symbol:
                        structure_image_url = f"{API_BASE_URL}/structure_image/{symbol}_{gene_id}_structure.png"
                        filtered_record['structure_image'] = f"![{symbol} Gene Structure]({structure_image_url})"
                filtered_result.append(filtered_record)

    return filtered_result if filtered_result else result


@mcp.tool()
async def fetch_gene_structure_details(entrezgene_id: str) -> dict:
    """
    Get detailed gene structure information (Transcripts, Exons, Introns)

    Use this for: transcript variants, exon/intron coordinates, splice sites, alternative splicing
    Use fetch_ncbi_info for: general gene info and structure visualization

    Args:
        entrezgene_id: NCBI Entrezgene ID (e.g., "839341"), supports semicolon-separated multiple IDs

    Returns:
        Dictionary with:
        - summary: Formatted markdown summary
        - result: List containing:
          - transcripts: List of transcript IDs (e.g., ["NM_001331246.1", ...])
          - exons: List of exon details with transcript, exon number, start, end, strand
          - introns: List of intron details with transcript, intron number, start, end, length
    """
    # Helper functions for parsing
    def parse_transcripts(text):
        if not text:
            return []
        return [t.replace('mRNA:', '').strip() for t in str(text).split(';') if t.strip()]

    def parse_exon_details(text):
        if not text:
            return []
        results = []
        for item in str(text).split('|'):
            item = item.strip()
            if not item:
                continue
            try:
                parts = item.split(':')
                if len(parts) >= 3:
                    transcript = parts[0]
                    exon_num = int(parts[1].replace('exon', ''))
                    pos_strand = parts[2]
                    pos, strand = pos_strand.rstrip(')').split('(')
                    start, end = pos.split('-')
                    results.append({
                        'transcript': transcript,
                        'exon': exon_num,
                        'start': int(start),
                        'end': int(end),
                        'strand': strand
                    })
            except (ValueError, IndexError):
                continue
        return results

    def parse_intron_details(text):
        if not text:
            return []
        results = []
        for item in str(text).split('|'):
            item = item.strip()
            if not item:
                continue
            try:
                parts = item.split(':')
                if len(parts) >= 3:
                    transcript = parts[0]
                    intron_num = int(parts[1].replace('intron', ''))
                    pos_len = parts[2]
                    pos, length_str = pos_len.split('(len=')
                    start, end = pos.split('-')
                    length = int(length_str.rstrip(')'))
                    results.append({
                        'transcript': transcript,
                        'intron': intron_num,
                        'start': int(start),
                        'end': int(end),
                        'length': length
                    })
            except (ValueError, IndexError):
                continue
        return results

    # Call HTTP API
    raw_result = await call_api("fetch_gene_structure_details", {"entrezgene_id": entrezgene_id})

    if not raw_result or (isinstance(raw_result, list) and len(raw_result) == 1 and "N/A" in raw_result[0]):
        return {"summary": f"No structure details found for {entrezgene_id}", "result": []}

    # Process each record
    processed_results = []
    summary_lines = [f"# Gene Structure Details for {entrezgene_id}\n"]

    for record in raw_result:
        if not isinstance(record, dict):
            continue

        eid = record.get('entrezgene_id', '')
        symbol = record.get('Symbol', '')

        transcripts = parse_transcripts(record.get('Transcripts'))
        exons = parse_exon_details(record.get('ExonDetails'))
        introns = parse_intron_details(record.get('IntronDetails'))

        processed_results.append({
            'entrezgene_id': eid,
            'symbol': symbol,
            'transcripts': transcripts,
            'exons': exons,
            'introns': introns
        })

        # Build summary
        summary_lines.append(f"\n## {symbol} ({eid})\n")
        summary_lines.append(f"### Transcripts ({len(transcripts)})")
        for t in transcripts[:10]:  # Limit display
            summary_lines.append(f"- {t}")
        if len(transcripts) > 10:
            summary_lines.append(f"- ... and {len(transcripts) - 10} more")

        # Group exons by transcript
        transcript_exons = {}
        for e in exons:
            t = e['transcript']
            if t not in transcript_exons:
                transcript_exons[t] = []
            transcript_exons[t].append(e)

        summary_lines.append(f"\n### Exons ({len(exons)} total across {len(transcript_exons)} transcripts)")
        for t, es in list(transcript_exons.items())[:3]:  # Show first 3 transcripts
            summary_lines.append(f"\n**{t}** ({len(es)} exons):")
            for e in es[:5]:  # Show first 5 exons
                summary_lines.append(f"- Exon {e['exon']}: {e['start']:,} - {e['end']:,} ({e['strand']})")
            if len(es) > 5:
                summary_lines.append(f"- ... and {len(es) - 5} more exons")

        # Group introns by transcript
        transcript_introns = {}
        for i in introns:
            t = i['transcript']
            if t not in transcript_introns:
                transcript_introns[t] = []
            transcript_introns[t].append(i)

        summary_lines.append(f"\n### Introns ({len(introns)} total)")
        for t, ins in list(transcript_introns.items())[:2]:  # Show first 2 transcripts
            summary_lines.append(f"\n**{t}** ({len(ins)} introns):")
            for i in ins[:3]:  # Show first 3 introns
                summary_lines.append(f"- Intron {i['intron']}: {i['start']:,} - {i['end']:,} (len={i['length']:,})")
            if len(ins) > 3:
                summary_lines.append(f"- ... and {len(ins) - 3} more introns")

    return {
        "summary": "\n".join(summary_lines),
        "result": processed_results
    }


@mcp.tool()
async def fetch_kegg_info(kegg_id: str) -> list:
    """
    Get KEGG functional annotation and pathway information

    Returns comprehensive information including pathways, enzyme functions, GO terms,
    orthology, disease associations, and external database links (27 fields total).

    Args:
        kegg_id: KEGG ID (e.g., "ath:AT1G01060"), supports semicolon-separated multiple IDs

    Returns:
        List of dictionaries with KEGG annotations including:
        PATHWAY_Names, PATHWAY_IDs, ENZYME, REACTION, GO terms, ORTHOLOGY,
        DISEASE, DEFINITION, DBLINKS, REFERENCE, etc.
    """
    return await call_api("fetch_kegg_info", {"kegg_id": kegg_id})


@mcp.tool()
async def fetch_uniprot_entry(gene_id: str, extentional_search: bool = False) -> str:
    """
    Get UniProt entry for gene IDs

    Args:
        gene_id: Gene ID to search, supports semicolon-separated multiple IDs
        extentional_search: Whether to perform extended search for unverified proteins

    Returns:
        Semicolon-separated uniprot_entry values
    """
    result = await call_api("fetch_uniprot_entry", {
        "gene_id": gene_id,
        "extended": str(extentional_search).lower()
    })
    return result.get("result", "N/A")


@mcp.tool()
async def fetch_uniprot_details(entries: str) -> list:
    """
    Get detailed UniProt information

    Args:
        entries: UniProt entry values, supports semicolon separation

    Returns:
        List of UniProt detail dictionaries with uniprot_link
    """
    result = await call_api("fetch_uniprot_details", {"entries": entries})
    if isinstance(result, list):
        for record in result:
            if isinstance(record, dict) and 'Entry' in record:
                record['uniprot_link'] = f"https://www.uniprot.org/uniprotkb/{record['Entry']}"
    return result


@mcp.tool()
async def fetch_interpro(interpro: str) -> list:
    """
    Get InterPro protein structure information

    Args:
        interpro: InterPro ID to search, supports semicolon-separated multiple IDs

    Returns:
        List of InterPro information dictionaries
    """
    return await call_api("fetch_interpro", {"interpro": interpro})


@mcp.tool()
async def fetch_pfam(pfam: str) -> list:
    """
    Get Pfam domain information

    Args:
        pfam: Pfam ID to search, supports semicolon-separated multiple IDs

    Returns:
        List of Pfam information dictionaries
    """
    return await call_api("fetch_pfam", {"pfam": pfam})


@mcp.tool()
async def fetch_pmid(query: str, max_results: int = 20) -> dict:
    """
    Search PubMed articles by query keywords

    Args:
        query: Query keywords for PubMed search
        max_results: Maximum number of results (default 20)

    Returns:
        Dictionary with article PMIDs and titles
    """
    if not ncbi_email or not ncbi_api_key:
        return {"error": "NCBI API not configured. Please set ncbi_config in config file or NCBI_EMAIL/NCBI_API_KEY environment variables."}

    try:
        async with httpx.AsyncClient(timeout=30.0) as client:
            # Search PubMed
            search_url = ncbi_base_url + "esearch.fcgi"
            params = {
                "db": "pubmed",
                "term": query,
                "retmax": max_results,
                "retmode": "xml",
                "email": ncbi_email,
                "api_key": ncbi_api_key
            }
            resp = await client.get(search_url, params=params)
            resp.raise_for_status()
            root = ET.fromstring(resp.text)
            pmids = [id_tag.text for id_tag in root.findall(".//Id")]

            if not pmids:
                return {"result": "No papers found"}

            # Fetch titles
            fetch_url = ncbi_base_url + "efetch.fcgi"
            params = {
                "db": "pubmed",
                "id": ",".join(pmids),
                "retmode": "xml",
                "email": ncbi_email,
                "api_key": ncbi_api_key
            }
            resp = await client.get(fetch_url, params=params)
            resp.raise_for_status()
            root = ET.fromstring(resp.text)

            papers = []
            for article in root.findall(".//PubmedArticle"):
                papers.append({
                    "pmid": article.findtext(".//PMID", default="N/A"),
                    "title": article.findtext(".//ArticleTitle", default="N/A")
                })
            return {"result": papers}
    except Exception as e:
        logger.error(f"PubMed search error: {e}")
        return {"error": str(e)}


@mcp.tool()
async def pubmed_advanced_search(
    query: str,
    page: int = 1,
    sort: str = "relevance",
    date_from: str = "",
    date_to: str = "",
    article_type: str = ""
) -> dict:
    """
    Advanced PubMed search with pagination (25 results per page)

    Supports PubMed advanced search syntax:
    - Field tags: [Title], [Author], [Journal], [MeSH], [Abstract], [Affiliation]
    - Boolean: AND, OR, NOT
    - Phrases: "exact phrase"
    - Wildcards: * (e.g., gene*)

    Examples:
    - "CRISPR[Title] AND plant[Title]" - Search CRISPR in plant titles
    - "Zhang Wei[Author] AND Arabidopsis" - Author + keyword
    - "drought tolerance[MeSH]" - MeSH term search
    - "Nature[Journal] AND 2024[Date - Publication]" - Journal + year

    Pagination:
    - page=1: First 25 results (default)
    - page=2: Results 26-50
    - page=3: Results 51-75
    - Use next_page value from response to get next page

    Args:
        query: PubMed search query (supports advanced syntax)
        page: Page number (default 1, each page has 25 results)
        sort: Sort order - "relevance" (default), "date", "author", "journal"
        date_from: Start date filter (YYYY/MM/DD or YYYY), empty string to skip
        date_to: End date filter (YYYY/MM/DD or YYYY), empty string to skip
        article_type: Filter by type - "review", "research", "clinical_trial", "meta_analysis", empty string to skip

    Returns:
        Dictionary with:
        - total_count: Total matching articles
        - current_page: Current page number
        - total_pages: Total number of pages available
        - results_per_page: Fixed at 25
        - returned_count: Number of results on this page
        - has_next_page: Whether more pages are available
        - next_page: Next page number (0 if no more pages)
        - has_prev_page: Whether previous page exists
        - prev_page: Previous page number (0 if on first page)
        - results: List of articles with PMID, Title, Authors, Journal, Year, Abstract snippet
        - query_translation: How PubMed interpreted your query
    """
    if not ncbi_email or not ncbi_api_key:
        return {"error": "NCBI API not configured"}

    # Fixed 25 results per page
    results_per_page = 25
    page = max(1, page)
    start = (page - 1) * results_per_page

    # Build search parameters
    search_params = {
        "db": "pubmed",
        "term": query,
        "retmax": results_per_page,
        "retstart": start,
        "retmode": "json",
        "email": ncbi_email,
        "api_key": ncbi_api_key,
        "usehistory": "y"
    }

    # Sort mapping
    sort_map = {
        "relevance": "relevance",
        "date": "pub_date",
        "author": "first_author",
        "journal": "journal"
    }
    if sort in sort_map:
        search_params["sort"] = sort_map[sort]

    # Date range filter
    if date_from:
        search_params["datetype"] = "pdat"
        search_params["mindate"] = date_from
    if date_to:
        search_params["datetype"] = "pdat"
        search_params["maxdate"] = date_to

    # Article type filter (append to query)
    type_filters = {
        "review": " AND review[pt]",
        "research": " AND journal article[pt]",
        "clinical_trial": " AND clinical trial[pt]",
        "meta_analysis": " AND meta-analysis[pt]"
    }
    if article_type and article_type in type_filters:
        search_params["term"] += type_filters[article_type]

    try:
        async with httpx.AsyncClient(timeout=60.0) as client:
            # Step 1: Search
            search_url = ncbi_base_url + "esearch.fcgi"
            resp = await client.get(search_url, params=search_params)
            resp.raise_for_status()
            search_result = resp.json()

            esearch = search_result.get("esearchresult", {})
            total_count = int(esearch.get("count", 0))
            pmids = esearch.get("idlist", [])
            query_translation = esearch.get("querytranslation", query)

            if not pmids:
                return {
                    "total_count": 0,
                    "current_page": page,
                    "total_pages": 0,
                    "results_per_page": results_per_page,
                    "returned_count": 0,
                    "has_next_page": False,
                    "next_page": 0,
                    "has_prev_page": False,
                    "prev_page": 0,
                    "query_translation": query_translation,
                    "results": [],
                    "message": "No articles found"
                }

            # Step 2: Fetch summaries (faster than full fetch)
            summary_url = ncbi_base_url + "esummary.fcgi"
            summary_params = {
                "db": "pubmed",
                "id": ",".join(pmids),
                "retmode": "json",
                "email": ncbi_email,
                "api_key": ncbi_api_key
            }
            resp = await client.get(summary_url, params=summary_params)
            resp.raise_for_status()
            summary_data = resp.json()

            # Parse results - NO ABSTRACT FETCHING (removed step 3)
            results = []
            uid_list = summary_data.get("result", {}).get("uids", [])
            for uid in uid_list:
                doc = summary_data.get("result", {}).get(uid, {})
                if not doc:
                    continue

                # Extract authors (first 3 + et al)
                authors = doc.get("authors", [])
                author_names = [a.get("name", "") for a in authors[:3]]
                if len(authors) > 3:
                    author_names.append("et al.")

                results.append({
                    "pmid": uid,
                    "title": doc.get("title", ""),
                    "authors": ", ".join(author_names),
                    "journal": doc.get("fulljournalname", doc.get("source", "")),
                    "year": doc.get("pubdate", "")[:4],
                    "pub_date": doc.get("pubdate", ""),
                    "doi": doc.get("elocationid", "").replace("doi: ", "") if "doi" in doc.get("elocationid", "").lower() else "",
                    "pubmed_link": f"https://pubmed.ncbi.nlm.nih.gov/{uid}/"
                })

            # Pagination info (no None values - MCP doesn't allow them)
            total_pages = (total_count + results_per_page - 1) // results_per_page if total_count > 0 else 1
            has_next_page = page < total_pages
            has_prev_page = page > 1

            return {
                "total_count": total_count,
                "current_page": page,
                "total_pages": total_pages,
                "results_per_page": results_per_page,
                "returned_count": len(results),
                "has_next_page": has_next_page,
                "next_page": page + 1 if has_next_page else 0,
                "has_prev_page": has_prev_page,
                "prev_page": page - 1 if has_prev_page else 0,
                "query_translation": query_translation,
                "results": results
            }

    except Exception as e:
        logger.error(f"PubMed advanced search error: {e}")
        return {"error": str(e)}    

@mcp.tool()
async def fetch_details(pmid: str) -> list:
    """
    Get detailed PubMed article information by PMID

    Args:
        pmid: PMIDs separated by commas

    Returns:
        List of detailed article information including title, authors, abstract, journal, year, DOI
    """
    if not ncbi_email or not ncbi_api_key:
        return [{"error": "NCBI API not configured. Please set ncbi_config in config file or NCBI_EMAIL/NCBI_API_KEY environment variables."}]

    if isinstance(pmid, str):
        pmid_list = [p.strip() for p in pmid.split(",") if p.strip()]
    else:
        pmid_list = pmid

    if not pmid_list:
        return [{"error": "No PMIDs provided"}]

    try:
        async with httpx.AsyncClient(timeout=30.0) as client:
            url = ncbi_base_url + "efetch.fcgi"
            params = {
                "db": "pubmed",
                "id": ",".join(pmid_list),
                "retmode": "xml",
                "email": ncbi_email,
                "api_key": ncbi_api_key
            }
            resp = await client.get(url, params=params)
            resp.raise_for_status()
            root = ET.fromstring(resp.text)

            results = []
            for article in root.findall(".//PubmedArticle"):
                info = {
                    "PMID": article.findtext(".//PMID"),
                    "Title": article.findtext(".//ArticleTitle"),
                    "Journal": article.findtext(".//Journal/Title"),
                    "Year": article.findtext(".//JournalIssue/PubDate/Year") or
                            article.findtext(".//JournalIssue/PubDate/MedlineDate"),
                }

                # Extract authors
                authors = []
                for author in article.findall(".//Author"):
                    ln = author.findtext("LastName")
                    fn = author.findtext("ForeName")
                    if ln and fn:
                        authors.append(f"{fn} {ln}")
                info["Authors"] = authors

                # Extract abstract
                abstract_texts = []
                for abst in article.findall(".//AbstractText"):
                    part = abst.text if abst.text else ""
                    label = abst.get("Label")
                    if label:
                        part = f"[{label}] {part}"
                    abstract_texts.append(part)
                info["Abstract"] = abstract_texts

                # Extract DOI
                for aid in article.findall(".//ArticleId"):
                    if aid.get("IdType") == "doi":
                        info["DOI"] = aid.text

                results.append(info)
            return results
    except Exception as e:
        logger.error(f"PubMed fetch error: {e}")
        return [{"error": str(e)}]


# ============================================================================
# STRING Database Tools (direct API calls)
# ============================================================================

if string_base_url:

    @mcp.tool(title="STRING: Resolves protein identifiers")
    async def string_resolve_proteins(
        proteins: Annotated[str, Field(description="Protein identifiers separated by %0d")],
        species: Annotated[str, Field(description="NCBI taxonomy ID")] = None,
        show_sequence: Annotated[str, Field(description="'1' to include sequences")] = None
    ) -> dict:
        """Maps protein identifiers to STRING metadata"""
        params = {"identifiers": proteins, "echo_query": 1, 'add_domains': 1}
        if species:
            params["species"] = species
        if show_sequence:
            params["add_sequence"] = show_sequence
        client = await get_string_client()
        log_call("/api/json/get_string_ids", params)
        results = await _post_json(client, "/api/json/get_string_ids", data=params)
        return {"mapped_proteins": results}

    @mcp.tool(title="STRING: Get interactions within query set")
    async def string_interactions_query_set(
        proteins: Annotated[str, Field(description="Protein identifiers separated by %0d")],
        species: Annotated[str, Field(description="NCBI taxonomy ID")] = None,
        required_score: Annotated[Optional[int], Field(description="Minimum confidence score 0-1000")] = None,
        network_type: Annotated[Optional[str], Field(description="'functional' or 'physical'")] = None,
        extend_network: Annotated[Optional[int], Field(description="Additional proteins to add")] = None,
    ) -> dict:
        """Retrieves interactions between query proteins"""
        params = {"identifiers": proteins}
        if species:
            params["species"] = species
        if required_score:
            params["required_score"] = required_score
        if network_type:
            params["network_type"] = network_type
        if extend_network:
            params["add_white_nodes"] = extend_network

        add_score_note = False
        if not required_score and len(proteins.lower().split("%d0")) <= 5:
            params["required_score"] = 0
            add_score_note = True

        client = await get_string_client()
        log_call("/api/json/network", params)
        results = await _post_json(client, "/api/json/network", data=params)
        if 'error' in results:
            return results
        results, add_size_note = truncate_network(results, required_score, 100)
        notes = []
        if add_score_note:
            notes.append("required_score lowered to 0.")
        if not results:
            notes.append("No interactions found.")
        if add_size_note:
            notes.append("List truncated to top 100.")
        return {"notes": notes, "network": results}

    @mcp.tool(title="STRING: Get all interaction partners")
    async def string_all_interaction_partners(
        identifiers: Annotated[str, Field(description="Protein identifiers separated by %0d")],
        species: Annotated[str, Field(description="NCBI taxonomy ID")] = None,
        required_score: Annotated[Optional[int], Field(description="Minimum score 0-1000")] = None,
        network_type: Annotated[Optional[str], Field(description="'functional' or 'physical'")] = None
    ) -> dict:
        """Retrieves all interaction partners for proteins"""
        params = {"identifiers": identifiers, "limit": 100}
        if species:
            params["species"] = species
        if required_score:
            params["required_score"] = required_score
        if network_type:
            params["network_type"] = network_type

        client = await get_string_client()
        log_call("/api/json/interaction_partners", params)
        results = await _post_json(client, "/api/json/interaction_partners", data=params)
        if 'error' in results:
            return results
        results_truncated, add_note = truncate_network(results, required_score, 100)
        notes = [f"Truncated to top 100."] if add_note else []
        return {"notes": notes, "interactions": results_truncated}

    @mcp.tool(title="STRING: Functional enrichment analysis")
    async def string_enrichment(
        proteins: Annotated[str, Field(description="Protein identifiers separated by %0d")],
        species: Annotated[Optional[str], Field(description="NCBI taxonomy ID")] = None
    ) -> dict:
        """Retrieves functional enrichment for proteins"""
        params = {"identifiers": proteins}
        if species:
            params["species"] = species

        client = await get_string_client()
        log_call("/api/json/enrichment", params)
        results = await _post_json(client, "/api/json/enrichment", data=params)
        if 'error' in results:
            return results
        results_truncated = truncate_enrichment(results, 'json')
        notes = []
        if not results_truncated:
            notes.append("No significant enrichment observed.")
        return {"notes": notes, "results": results_truncated}

    @mcp.tool(title="STRING: Get functional annotations")
    async def string_functional_annotation(
        identifiers: Annotated[str, Field(description="Protein identifiers separated by %0d")],
        species: Annotated[str, Field(description="NCBI taxonomy ID")] = None,
    ) -> dict:
        """Retrieves functional annotations for proteins"""
        params = {"identifiers": identifiers}
        if species:
            params["species"] = species

        client = await get_string_client()
        log_call("/api/json/functional_annotation", params)
        results = await _post_json(client, "/api/json/functional_annotation", data=params)
        if 'error' in results:
            return results
        # Sort and truncate
        if isinstance(results, list):
            results = sorted(results, key=lambda x: x.get("ratio_in_set", 0), reverse=True)[:50]
        return {"results": results}

    @mcp.tool(title="STRING: Query species")
    async def string_query_species(
        species_text: Annotated[str, Field(description="Species name or taxon ID to search")]
    ) -> dict:
        """Search for species in STRING"""
        params = {"species_text": species_text, 'limit': 50, 'add_sps': 't'}
        client = await get_string_client()
        log_call("/api/json/query_species_names", params)
        results = await _post_json(client, "/api/json/query_species_names", data=params)
        return {"results": results}

    @mcp.tool(title="STRING: Get network image URL")
    async def string_visual_network(
        proteins: Annotated[str, Field(description="Protein identifiers")],
        species: Annotated[str, Field(description="NCBI taxonomy ID")] = None,
        extend_network: Annotated[Optional[int], Field(description="Add nodes")] = None,
        required_score: Annotated[Optional[int], Field(description="Score threshold 0-1000")] = None,
        network_type: Annotated[Optional[str], Field(description="'functional' or 'physical'")] = None,
        network_flavor: Annotated[Optional[str], Field(description="'evidence', 'confidence', 'actions'")] = None,
        hide_disconnected_nodes: Annotated[Optional[int], Field(description="1 to hide")] = None,
    ) -> dict:
        """Gets STRING network image URL"""
        params = {"identifiers": proteins}
        if species:
            params["species"] = species
        if extend_network:
            params["add_white_nodes"] = extend_network
        if required_score:
            params["required_score"] = required_score
        if network_type:
            params["network_type"] = network_type
        if network_flavor:
            params["network_flavor"] = network_flavor
        if hide_disconnected_nodes:
            params["hide_disconnected_nodes"] = hide_disconnected_nodes

        if not required_score and len(proteins.lower().split("%d0")) <= 5:
            params["required_score"] = 0

        client = await get_string_client()
        log_call("/api/json/network_image_url", params)
        results = await _post_json(client, "/api/json/network_image_url", data=params)
        return {"notes": ["Embed as markdown image."], "image_url": results}

    @mcp.tool(title="STRING: Get interactive network link")
    async def string_network_link(
        proteins: Annotated[str, Field(description="Protein identifiers")],
        species: Annotated[str, Field(description="NCBI taxonomy ID")] = None,
        extend_network: Annotated[Optional[int], Field(description="Add nodes")] = None,
        required_score: Annotated[Optional[int], Field(description="Score threshold")] = None,
        network_flavor: Annotated[Optional[str], Field(description="Edge style")] = None,
        network_type: Annotated[Optional[str], Field(description="'functional' or 'physical'")] = None,
        hide_disconnected_nodes: Annotated[Optional[int], Field(description="1 to hide")] = None,
    ) -> dict:
        """Gets interactive STRING network link"""
        params = {"identifiers": proteins}
        if species:
            params["species"] = species
        if extend_network:
            params["add_white_nodes"] = extend_network
        if required_score:
            params["required_score"] = required_score
        if network_flavor:
            params["network_flavor"] = network_flavor
        if network_type:
            params["network_type"] = network_type
        if hide_disconnected_nodes:
            params["hide_disconnected_nodes"] = hide_disconnected_nodes

        if not required_score and len(proteins.lower().split("%d0")) <= 5:
            params["required_score"] = 0

        client = await get_string_client()
        log_call("/api/json/get_link", params)
        results = await _post_json(client, "/api/json/get_link", data=params)
        return {"notes": ["Embed as markdown link."], "results": results}

    @mcp.tool(title="STRING: Get homologs")
    async def string_homology(
        proteins: Annotated[str, Field(description="Protein identifiers separated by %0d")],
        species: Annotated[str, Field(description="NCBI taxonomy ID")] = None,
        species_b: Annotated[Optional[str], Field(description="Target species IDs")] = None
    ) -> dict:
        """Retrieves protein homologs"""
        params = {"identifiers": proteins}
        if species:
            params["species"] = species
        if species_b:
            params["species_b"] = species_b

        client = await get_string_client()
        log_call("/api/json/homology_all", params)
        results = await _post_json(client, "/api/json/homology_all", data=params)
        return {"results": results}

    @mcp.tool(title="STRING: Get interaction evidence links")
    async def string_interaction_evidence(
        identifier_a: Annotated[str, Field(description="Protein A identifier")],
        identifiers_b: Annotated[str, Field(description="Protein B identifiers separated by %0d")],
        species: Annotated[str, Field(description="NCBI taxonomy ID")] = None,
    ) -> dict:
        """Gets links to interaction evidence pages"""
        identifiers_b = identifiers_b.replace('%0D', '%0d')
        output = [f"{string_base_url}/interaction/{identifier_a}/{b}?species={species}&suppress_disambiguation=1"
                  for b in identifiers_b.split("%0d")]
        return {"notes": ["Embed as markdown links."], "results": output}

    @mcp.tool(title="STRING: Get enrichment image URL")
    async def string_enrichment_image_url(
        identifiers: Annotated[str, Field(description="Protein identifiers separated by %0d")],
        species: Annotated[str, Field(description="NCBI taxonomy ID")] = None,
        category: Annotated[Optional[str], Field(description="Term category")] = None,
    ) -> dict:
        """Gets enrichment figure URL"""
        params = {"identifiers": identifiers}
        if species:
            params["species"] = species
        if category:
            params["category"] = category

        client = await get_string_client()
        log_call("/api/json/enrichment_image_url", params)
        results = await _post_json(client, "/api/json/enrichment_image_url", data=params)
        return {"notes": ["Embed as markdown image."], "results": results}

    @mcp.tool(title="STRING: PPI enrichment")
    async def string_ppi_enrichment(
        identifiers: Annotated[str, Field(description="Protein identifiers separated by %0d")],
        species: Annotated[str, Field(description="NCBI taxonomy ID")] = None,
        required_score: Annotated[Optional[int], Field(description="Minimum score")] = None,
    ) -> dict:
        """Tests if network is enriched in protein-protein interactions"""
        params = {"identifiers": identifiers}
        if species:
            params["species"] = species
        if required_score:
            params["required_score"] = required_score

        client = await get_string_client()
        log_call("/api/json/ppi_enrichment", params)
        results = await _post_json(client, "/api/json/ppi_enrichment", data=params)
        return {"results": results}

    @mcp.tool(title="STRING: Get proteins for functional term")
    async def string_proteins_for_term(
        term_text: Annotated[str, Field(description="Functional term or text")],
        species: Annotated[str, Field(description="NCBI taxonomy ID")] = "9606"
    ) -> dict:
        """Retrieves proteins annotated with a functional term"""
        params = {"term_text": term_text, "species": species}
        client = await get_string_client()
        log_call("/api/json/functional_terms", params)
        results = await _post_json(client, "/api/json/functional_terms", data=params)
        if 'error' in results:
            return results
        # Truncate
        if isinstance(results, list):
            for i, row in enumerate(results[:10]):
                cutoff = 100 if i < 3 else 25
                if len(row.get('preferredNames', [])) > cutoff:
                    row['preferredNames'] = row['preferredNames'][:cutoff] + ["..."]
                    row['stringIds'] = row['stringIds'][:cutoff] + ["..."]
            results = results[:10]
        return {"results": results}

    @mcp.tool(title="STRING: Sequence search")
    async def string_sequence_search(
        sequences: Annotated[str, Field(description="Protein sequences in FASTA format")],
        species: Annotated[str, Field(description="NCBI taxonomy ID")] = "9606"
    ) -> dict:
        """Searches STRING by amino acid sequence"""
        params = {"sequences": sequences, "species": species}
        # Use longer timeout for sequence search
        client = await get_string_client()
        log_call("/api/json/similarity_search", params)
        results = await _post_json(client, "/api/json/similarity_search", data=params)
        # Truncate to 50
        if isinstance(results, list) and len(results) > 50:
            results = sorted(results, key=lambda x: x.get("bitscore", 0), reverse=True)[:50]
        return {"results": results}

    @mcp.tool(title="STRING: Help/FAQ")
    async def string_help(
        topic: Annotated[str, Field(description="Help topic")]
    ) -> dict:
        """Provides help text for STRING features"""
        key = topic.lower()
        if key not in HELP_TOPICS:
            return {"error": f"Unknown topic. Available: {', '.join(HELP_TOPICS.keys())}"}
        return {"topic": key, "text": HELP_TOPICS[key]}


# ============================================================================
# Service State
# ============================================================================

class ServiceState:
    def __init__(self):
        self.initialized = False

    def set_initialized(self):
        self.initialized = True

    def is_initialized(self):
        return self.initialized


service_state = ServiceState()


@mcp.tool()
async def check_initialization():
    """Check if service is initialized"""
    if not service_state.is_initialized():
        raise Exception("Service not initialized.")
    return {"status": "initialized"}


def check_keyboard_input():
    """Thread for keyboard input"""
    while True:
        try:
            user_input = input()
            if user_input.lower() in ['q', 'quit', 'exit']:
                logger.info("Exiting...")
                os._exit(0)
        except:
            pass
        time.sleep(0.1)


# ============================================================================
# Main Entry Point
# ============================================================================

if __name__ == "__main__":
    keyboard_thread = threading.Thread(target=check_keyboard_input, daemon=True)
    keyboard_thread.start()

    logger.info("=" * 60)
    logger.info("Plant Gene & STRING Database MCP Server")
    logger.info(f"API Backend: {API_BASE_URL}")
    if string_base_url:
        logger.info(f"STRING API: {string_base_url}")
    logger.info("=" * 60)

    # Initialize locks before starting server (prevents race condition)
    _ensure_locks_initialized()

    time.sleep(0.5)
    service_state.set_initialized()
    logger.info("Service initialized.")

    logger.info("Transport: HTTP")
    logger.info("Host: 0.0.0.0")
    logger.info("Port: 1154")
    logger.info("Enter 'q' to exit.")
    logger.info("=" * 60)

    mcp.run(transport="http", port=1154, host="0.0.0.0")
