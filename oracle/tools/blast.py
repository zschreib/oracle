"""
BLAST tool for ORACLE.

Hits the NCBI BLASTp REST API and parses the top hit into an
AnnotationEvidence object. Uses the two-step NCBI BLAST API pattern:
  1. Submit job → get RID (request ID)
  2. Poll until complete → parse results

NCBI rate limits unauthenticated requests to ~3/second. We add
conservative sleep intervals to stay well under that. If you have
an NCBI API key, set it in the environment as NCBI_API_KEY and
the rate limit rises to 10/second.

Coverage is computed as (alignment length / query length) and is
the most important filter alongside e-value. A high-scoring hit
with low coverage (e.g. only a domain matched) is weaker evidence
than the raw e-value suggests.
"""

import os
import time
import requests
from xml.etree import ElementTree as ET

from .base import AnnotationTool
from ..models import AnnotationEvidence, ConfidenceTier


# NCBI BLAST REST endpoints
_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"

# Confidence thresholds based on e-value and coverage.
# These are conservative: better to call MODERATE and escalate
# than to call HIGH and stop too early.
_HIGH_EVALUE    = 1e-10  # strong homology
_MODERATE_EVALUE = 1e-3  # weak but meaningful signal
_HIGH_COVERAGE  = 0.70   # at least 70% of query covered by hit
_MIN_COVERAGE   = 0.40   # below this, hit is probably a partial domain match


class BlastTool(AnnotationTool):
    """
    BLASTp search against NCBI nr database via REST API.

    The agent should call this first for any unannotated protein sequence.
    Results with e-value > 1e-3 or coverage < 0.40 are returned as LOW
    confidence and should trigger HMMer escalation.
    """

    def __init__(self, database: str = "nr", max_hits: int = 5,
                 poll_interval: int = 15, max_wait: int = 300):
        """
        Args:
            database:      NCBI database to search. 'nr' is the default
                           non-redundant protein database. 'refseq_protein'
                           is faster and cleaner but less comprehensive.
            max_hits:      Number of top hits to retrieve. We use the top
                           hit for annotation but log others for context.
            poll_interval: Seconds between status checks while job runs.
            max_wait:      Maximum seconds to wait before timing out.
        """
        self.database = database
        self.max_hits = max_hits
        self.poll_interval = poll_interval
        self.max_wait = max_wait
        self.api_key = os.environ.get("NCBI_API_KEY", "")

    @property
    def name(self) -> str:
        return "blast_nr"

    @property
    def tool_schema(self) -> dict:
        return {
            "name": self.name,
            "description": (
                "Run BLASTp homology search against the NCBI nr database. "
                "Use this as the first step for any unannotated protein sequence. "
                "Returns the top hit with e-value, percent identity, and coverage. "
                "Do NOT use if you already have a strong BLAST result (e-value < 1e-10, "
                "coverage > 0.70) from a previous step. "
                "Takes 30-120 seconds due to NCBI queue times."
            ),
            "input_schema": {
                "type": "object",
                "properties": {
                    "sequence": {
                        "type": "string",
                        "description": "Amino acid sequence string, no FASTA header, no whitespace."
                    },
                    "database": {
                        "type": "string",
                        "description": "NCBI database. Default 'nr'. Use 'refseq_protein' for faster results.",
                        "default": "nr"
                    }
                },
                "required": ["sequence"]
            }
        }

    def run(self, sequence: str, database: str = None, **kwargs) -> AnnotationEvidence:
        """
        Submit BLASTp job, poll for completion, parse top hit.

        Never raises: any failure returns AnnotationEvidence with
        ConfidenceTier.UNKNOWN and the error message in reasoning.
        """
        db = database or self.database
        sequence = sequence.strip().replace("\n", "").replace(" ", "")

        try:
            rid = self._submit(sequence, db)
        except Exception as e:
            return self._failure(f"Failed to submit BLAST job: {e}")

        try:
            raw_xml = self._poll(rid)
        except Exception as e:
            return self._failure(f"BLAST polling failed or timed out: {e}")

        try:
            return self._parse(raw_xml, sequence)
        except Exception as e:
            return self._failure(f"Failed to parse BLAST results: {e}")

    def _submit(self, sequence: str, database: str) -> str:
        """Submit BLASTp job and return the RID for polling."""
        params = {
            "CMD": "Put",
            "PROGRAM": "blastp",
            "DATABASE": database,
            "QUERY": sequence,
            "FORMAT_TYPE": "XML",
            "HITLIST_SIZE": self.max_hits,
            # EMAIL is required by NCBI for automated requests
            "EMAIL": "oracle-annotation-agent@placeholder.com",
        }
        if self.api_key:
            params["api_key"] = self.api_key

        response = requests.post(_BLAST_URL, data=params, timeout=30)
        response.raise_for_status()

        # RID is embedded in the HTML response as "RID = <value>"
        for line in response.text.splitlines():
            if line.strip().startswith("RID ="):
                return line.strip().split("=")[1].strip()

        raise ValueError("RID not found in BLAST submission response")

    def _poll(self, rid: str) -> str:
        """Poll NCBI until job completes, then return raw XML results."""
        elapsed = 0
        while elapsed < self.max_wait:
            time.sleep(self.poll_interval)
            elapsed += self.poll_interval

            response = requests.get(_BLAST_URL, params={
                "CMD": "Get",
                "RID": rid,
                "FORMAT_TYPE": "XML",
            }, timeout=30)
            response.raise_for_status()

            # NCBI signals completion by the absence of "Status=WAITING"
            # and presence of BlastOutput XML
            if "Status=WAITING" in response.text:
                continue
            if "Status=FAILED" in response.text:
                raise RuntimeError(f"BLAST job {rid} failed on NCBI servers")
            if "<BlastOutput>" in response.text:
                return response.text

        raise TimeoutError(
            f"BLAST job {rid} did not complete within {self.max_wait}s"
        )

    def _parse(self, xml_text: str, query_sequence: str) -> AnnotationEvidence:
        """
        Parse BLAST XML output and return AnnotationEvidence.

        Extracts the top hit only. Coverage is computed from alignment
        length relative to query length, not subject length, because we
        care how much of our unknown sequence is explained by the hit.
        """
        root = ET.fromstring(xml_text)
        query_length = len(query_sequence)

        # Navigate BLAST XML structure to first hit
        iterations = root.findall(".//Iteration")
        if not iterations:
            return self._no_hits()

        hits = iterations[0].findall(".//Hit")
        if not hits:
            return self._no_hits()

        top_hit = hits[0]
        hsp = top_hit.find(".//Hsp")  # high scoring pair (alignment block)
        if hsp is None:
            return self._no_hits()

        # Extract alignment metrics
        evalue    = float(hsp.findtext("Hsp_evalue", default="1.0"))
        bitscore  = float(hsp.findtext("Hsp_bit-score", default="0.0"))
        align_len = int(hsp.findtext("Hsp_align-len", default="0"))
        identity  = int(hsp.findtext("Hsp_identity", default="0"))
        coverage  = align_len / query_length if query_length > 0 else 0.0
        pct_id    = identity / align_len * 100 if align_len > 0 else 0.0

        hit_id   = top_hit.findtext("Hit_id", default="unknown")
        hit_def  = top_hit.findtext("Hit_def", default="hypothetical protein")

        # Assign confidence tier based on e-value and coverage together.
        # Coverage matters: a 1e-50 hit covering 20% of the query is still
        # weak evidence for the full protein's function.
        confidence, reasoning = self._assign_confidence(
            evalue, coverage, pct_id, hit_def
        )

        # Clean up the annotation label from NCBI's verbose hit description
        annotation = self._clean_annotation(hit_def)

        return AnnotationEvidence(
            tool_name=self.name,
            annotation=annotation,
            confidence=confidence,
            score=evalue,
            coverage=coverage,
            hit_id=hit_id,
            hit_description=hit_def,
            reasoning=reasoning,
            raw_output=xml_text[:500]  # truncate for storage
        )

    def _assign_confidence(self, evalue: float, coverage: float,
                           pct_id: float, hit_def: str) -> tuple:
        """
        Map e-value + coverage to a ConfidenceTier with reasoning text.

        Returns (ConfidenceTier, reasoning_string).
        """
        # Hypothetical protein hits are weak even at good e-values.
        # Only flag as hypothetical if the hit description contains these
        # terms AND does not also contain a real functional keyword.
        # This prevents flagging "3'-5' exonuclease" as hypothetical
        # just because the NCBI reference entry has generic text.
        weak_terms = ["hypothetical", "uncharacterized", "unknown function", "putative"]
        functional_terms = ["polymerase", "exonuclease", "kinase", "protease",
                           "helicase", "lyase", "synthase", "transferase",
                           "reductase", "ligase", "hydrolase", "oxidase",
                           "fiber", "capsid", "tail", "integrase", "recombinase"]
        hit_lower = hit_def.lower()
        is_hypothetical = (
            any(term in hit_lower for term in weak_terms)
            and not any(term in hit_lower for term in functional_terms)
        )

        if evalue <= _HIGH_EVALUE and coverage >= _HIGH_COVERAGE and not is_hypothetical:
            tier = ConfidenceTier.HIGH
            reason = (
                f"Strong homology: e-value {evalue:.2e}, {coverage:.0%} query coverage, "
                f"{pct_id:.1f}% identity. Hit is a named protein. "
                f"Annotation is well-supported without further tool escalation."
            )
        elif evalue == 0.0 and coverage >= _HIGH_COVERAGE:
            # e-value of exactly 0.0 means the match is so perfect NCBI
            # cannot compute a meaningful probability. Always HIGH regardless
            # of hit description since the sequence identity is essentially 100%.
            tier = ConfidenceTier.HIGH
            reason = (
                f"Perfect match: e-value 0.0, {coverage:.0%} coverage, "
                f"{pct_id:.1f}% identity. Sequence is identical or near-identical "
                f"to a database entry. Annotation is definitive."
            )
        elif evalue <= _MODERATE_EVALUE and coverage >= _MIN_COVERAGE:
            tier = ConfidenceTier.MODERATE
            if is_hypothetical:
                reason = (
                    f"Hit is a hypothetical/uncharacterized protein (e-value {evalue:.2e}, "
                    f"{coverage:.0%} coverage). Meaningful signal but annotation is "
                    f"speculative. HMMer domain search recommended."
                )
            else:
                reason = (
                    f"Moderate homology: e-value {evalue:.2e}, {coverage:.0%} coverage, "
                    f"{pct_id:.1f}% identity. Hit has a functional name but coverage or "
                    f"e-value is below high-confidence threshold. HMMer recommended."
                )
        else:
            tier = ConfidenceTier.LOW
            reason = (
                f"Weak or no meaningful homology: e-value {evalue:.2e}, "
                f"{coverage:.0%} coverage. BLAST evidence is insufficient for "
                f"functional annotation. Escalate to HMMer domain search."
            )

        return tier, reason

    def _clean_annotation(self, hit_def: str) -> str:
        """
        Extract a clean annotation label from NCBI's hit description.

        NCBI descriptions are verbose: "tail fiber protein [Escherichia phage T4] >..."
        We want just the functional part before the organism bracket.
        """
        # Take only the first hit description if multiple are concatenated
        label = hit_def.split(">")[0].strip()
        # Remove organism name in brackets
        if "[" in label:
            label = label[:label.index("[")].strip()
        # Truncate if still very long
        return label[:80] if len(label) > 80 else label

    def _no_hits(self) -> AnnotationEvidence:
        """Return UNKNOWN evidence when BLAST finds no hits."""
        return AnnotationEvidence(
            tool_name=self.name,
            annotation="No significant hits",
            confidence=ConfidenceTier.UNKNOWN,
            reasoning=(
                "BLASTp returned no hits above threshold. Sequence has no "
                "detectable homology in NCBI nr. Escalate to HMMer domain "
                "search to check for conserved domain signatures."
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