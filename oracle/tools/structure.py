"""
RCSB PDB sequence similarity search tool for ORACLE.

Uses the RCSB PDB Search API (MMseqs2 backend) to find structurally
and functionally characterized proteins similar to the query sequence.
This is the third tier in the ORACLE evidence ladder, called when both
BLAST against nr and InterProScan domain search fail to produce
high-confidence functional annotations.

The key difference from BLAST (tier 1) is the database and sensitivity:
- BLAST searches nr (all sequences, many uncharacterized)
- This tool searches PDB (experimentally determined structures only,
  all entries have at minimum a structural characterization)

A hit against PDB with known structure provides stronger functional
inference than a hit against nr hypothetical proteins, because PDB
entries are experimentally validated and structurally characterized.

API: RCSB PDB Search API v2
Docs: https://search.rcsb.org/
No authentication required. Synchronous response, no polling needed.
"""

import requests

from .base import AnnotationTool
from ..models import AnnotationEvidence, ConfidenceTier


_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
_ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity"

_HIGH_IDENTITY = 0.70
_MOD_IDENTITY  = 0.30


class RCSBStructureTool(AnnotationTool):
    """
    RCSB PDB sequence similarity search via MMseqs2.

    Searches experimentally characterized PDB structures for sequence
    similarity. Called when BLAST and InterProScan both fail to produce
    high-confidence functional annotations.
    """

    def __init__(self, identity_cutoff: float = 0.20, max_hits: int = 5):
        self.identity_cutoff = identity_cutoff
        self.max_hits = max_hits

    @property
    def name(self) -> str:
        return "rcsb_pdb"

    @property
    def tool_schema(self) -> dict:
        return {
            "name": self.name,
            "description": (
                "Search the RCSB Protein Data Bank for structurally and "
                "functionally similar proteins using sequence similarity. "
                "Use this as the last resort when both blast_nr and hmmer_pfam "
                "return LOW or UNKNOWN confidence. "
                "Unlike BLAST against nr, this searches only experimentally "
                "determined structures with validated functional annotations. "
                "A hit here provides stronger functional inference than a "
                "hypothetical protein match from BLAST. "
                "Takes 5-15 seconds."
            ),
            "input_schema": {
                "type": "object",
                "properties": {
                    "sequence": {
                        "type": "string",
                        "description": (
                            "Amino acid sequence string, no FASTA header, "
                            "no whitespace."
                        )
                    }
                },
                "required": ["sequence"]
            }
        }

    def run(self, sequence: str, **kwargs) -> AnnotationEvidence:
        """
        Search RCSB PDB and return the best annotated hit.

        First searches experimental structures only. If no hits found,
        falls back to AlphaFold Computed Structure Models (CSMs) which
        cover 214 million proteins including environmental sequences.
        """
        sequence = sequence.strip().replace("\n", "").replace(" ", "")

        try:
            hits = self._search(sequence, content_type="experimental")
        except Exception as e:
            return self._failure(f"RCSB PDB search failed: {e}")

        source = "PDB (experimental)"

        # Fallback to AlphaFold CSMs if no experimental hits
        if not hits:
            try:
                hits = self._search(sequence, content_type="computational",
                                    timeout=60)
                source = "AlphaFold CSM"
            except Exception as e:
                return self._failure(f"RCSB CSM fallback search failed: {e}")

        if not hits:
            return self._no_hits()

        try:
            return self._annotate_top_hit(hits, sequence, source)
        except Exception as e:
            return self._failure(f"Failed to annotate RCSB hit: {e}")

    def _search(self, sequence: str,
                content_type: str = "experimental",
                timeout: int = 30) -> list:
        """
        Submit sequence to RCSB Search API and return result set.

        Args:
            sequence:     Amino acid sequence string
            content_type: 'experimental' for PDB structures only,
                          'computational' for AlphaFold CSMs only
        """
        query = {
            "query": {
                "type": "terminal",
                "service": "sequence",
                "parameters": {
                    "evalue_cutoff": 1,
                    "identity_cutoff": self.identity_cutoff,
                    "sequence_type": "protein",
                    "value": sequence
                }
            },
            "request_options": {
                "scoring_strategy": "sequence",
                "results_content_type": [content_type],
                "paginate": {"start": 0, "rows": self.max_hits}
            },
            "return_type": "polymer_entity"
        }

        response = requests.post(
            _SEARCH_URL,
            json=query,
            headers={"Content-Type": "application/json"},
            timeout=timeout
        )

        # RCSB returns 204 No Content when there are no hits
        if response.status_code == 204 or not response.text.strip():
            return []

        response.raise_for_status()
        return response.json().get("result_set", [])

    def _annotate_top_hit(self, hits: list, query_sequence: str,
                          source: str = "PDB (experimental)") -> AnnotationEvidence:
        """Fetch functional annotation for the top hit and build evidence."""
        top_hit = hits[0]
        identifier = top_hit.get("identifier", "")
        score = float(top_hit.get("score", 0.0))

        if "_" in identifier:
            pdb_id, entity_id = identifier.split("_", 1)
        else:
            pdb_id, entity_id = identifier, "1"

        entity_info = self._fetch_entity_info(pdb_id, entity_id)

        description = "Unknown protein"
        organism = ""

        if entity_info:
            poly_entity = entity_info.get("rcsb_polymer_entity", {})
            description = poly_entity.get("pdbx_description", description)
            src_org = entity_info.get("rcsb_entity_source_organism", [{}])
            if src_org:
                organism = src_org[0].get("ncbi_scientific_name", "")

        annotation = f"{description} [{organism}]" if organism else description
        annotation = annotation[:100] if len(annotation) > 100 else annotation

        confidence, reasoning = self._assign_confidence(
            score, description, pdb_id, source
        )

        return AnnotationEvidence(
            tool_name=self.name,
            annotation=annotation,
            confidence=confidence,
            score=score,
            coverage=score,
            hit_id=identifier,
            hit_description=f"{description} | {pdb_id} | {organism} | {source}",
            reasoning=reasoning,
            raw_output=str(hits[:2])
        )

    def _fetch_entity_info(self, pdb_id: str, entity_id: str) -> dict:
        """Fetch polymer entity metadata from RCSB Data API."""
        try:
            response = requests.get(
                f"{_ENTITY_URL}/{pdb_id}/{entity_id}",
                timeout=15
            )
            if response.status_code == 200:
                return response.json()
        except Exception:
            pass
        return {}

    def _assign_confidence(self, score: float, description: str,
                           pdb_id: str, source: str = "PDB") -> tuple:
        """Map RCSB similarity score to ConfidenceTier with reasoning."""
        is_vague = any(
            term in description.lower()
            for term in ["hypothetical", "uncharacterized", "unknown",
                         "putative", "unnamed"]
        )

        if score >= _HIGH_IDENTITY and not is_vague:
            tier = ConfidenceTier.HIGH
            reason = (
                f"Strong {source} hit: {score:.0%} sequence identity to {pdb_id}. "
                f"Experimentally characterized structure with validated function. "
                f"Annotation is well-supported."
            )
        elif score >= _MOD_IDENTITY and not is_vague:
            tier = ConfidenceTier.MODERATE
            reason = (
                f"Moderate {source} hit: {score:.0%} sequence identity to {pdb_id}. "
                f"Likely shares structural fold and general function. "
                f"Experimental validation recommended for specific activity."
            )
        else:
            tier = ConfidenceTier.LOW
            reason = (
                f"Weak {source} hit: {score:.0%} identity to {pdb_id}. "
                f"Below threshold for reliable functional inference. "
                f"Sequence may be genuine dark matter."
            )

        return tier, reason

        if score >= _HIGH_IDENTITY and not is_vague:
            tier = ConfidenceTier.HIGH
            reason = (
                f"Strong PDB hit: {score:.0%} sequence identity to {pdb_id}. "
                f"Experimentally characterized structure with validated function. "
                f"Annotation is well-supported."
            )
        elif score >= _MOD_IDENTITY and not is_vague:
            tier = ConfidenceTier.MODERATE
            reason = (
                f"Moderate PDB hit: {score:.0%} sequence identity to {pdb_id}. "
                f"Likely shares structural fold and general function. "
                f"Experimental validation recommended for specific activity."
            )
        else:
            tier = ConfidenceTier.LOW
            reason = (
                f"Weak PDB hit: {score:.0%} identity to {pdb_id}. "
                f"Below threshold for reliable functional inference. "
                f"Sequence may be genuine dark matter."
            )

        return tier, reason

    def _no_hits(self) -> AnnotationEvidence:
        """Return UNKNOWN evidence when RCSB finds no hits."""
        return AnnotationEvidence(
            tool_name=self.name,
            annotation="No PDB hits",
            confidence=ConfidenceTier.UNKNOWN,
            reasoning=(
                "RCSB PDB sequence search found no hits above threshold. "
                "Sequence has no detectable similarity to any experimentally "
                "characterized protein structure. This is genuine dark matter: "
                "no homologs in sequence databases, no conserved domains, "
                "and no structurally characterized relatives."
            )
        )

    def _failure(self, message: str) -> AnnotationEvidence:
        """Return UNKNOWN evidence on any tool failure."""
        return AnnotationEvidence(
            tool_name=self.name,
            annotation="Tool failure",
            confidence=ConfidenceTier.UNKNOWN,
            reasoning=message
        )