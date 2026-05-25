"""
Domain annotation tool for ORACLE using the EBI InterProScan REST API.

InterProScan searches multiple profile HMM and pattern databases simultaneously
including Pfam, TIGRFAM, Gene3D, PANTHER, CDD, SUPERFAMILY, and PRINTS.
This gives broader domain coverage than HMMer against Pfam alone and is more
reliable since InterProScan is a primary EBI service with SLA support.

This is the second tier in the ORACLE evidence ladder, called when BLAST
returns a hypothetical protein hit, low coverage, or weak e-value.

API: EBI Job Dispatcher REST for InterProScan 5
Docs: https://www.ebi.ac.uk/jdispatcher/pfa/iprscan5

Three-step pattern:
  1. POST form data to /run -> returns plain text job ID
  2. GET /status/{jobId} -> poll until FINISHED
  3. GET /result/{jobId}/json -> fetch structured results

Confidence assignment:
- InterProScan returns hits from multiple databases with e-values and scores
- We rank by database reliability: Pfam > TIGRFAM > Gene3D > others
- A hit with a curated description above threshold is HIGH confidence
- A DUF or uncharacterized domain is MODERATE
- No hits above threshold is UNKNOWN, escalate to RCSB PDB
"""

import time
import requests

from .base import AnnotationTool
from ..models import AnnotationEvidence, ConfidenceTier


_BASE       = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
_RUN_URL    = f"{_BASE}/run"
_STATUS_URL = f"{_BASE}/status"
_RESULT_URL = f"{_BASE}/result"

# Database priority order for selecting the best hit when multiple databases
# return results. More curated databases rank higher.
_DB_PRIORITY = {
    "Pfam": 10, "TIGRFAM": 9, "Gene3D": 8,
    "PANTHER": 7, "CDD": 6, "SUPERFAMILY": 5,
    "PRINTS": 4, "ProSiteProfiles": 3, "HAMAP": 3,
}

# Confidence thresholds
_HIGH_EVALUE   = 1e-5
_MOD_EVALUE    = 1e-2
_HIGH_COVERAGE = 0.60
_MIN_COVERAGE  = 0.25


class HmmerTool(AnnotationTool):
    """
    InterProScan domain search via EBI Job Dispatcher REST API.

    Despite the class name kept as HmmerTool for registry compatibility,
    this uses InterProScan which runs multiple domain databases including
    Pfam. Called when BLAST returns weak or hypothetical hits.
    """

    def __init__(self, poll_interval: int = 10, max_wait: int = 300,
                 email: str = "oracle-agent@example.com"):
        """
        Args:
            poll_interval: Seconds between status checks. InterProScan
                           typically takes 30-90 seconds.
            max_wait:      Maximum seconds to wait before timing out.
            email:         Required by EBI for job submission tracking.
        """
        self.poll_interval = poll_interval
        self.max_wait = max_wait
        self.email = email

    @property
    def name(self) -> str:
        return "hmmer_pfam"

    @property
    def tool_schema(self) -> dict:
        return {
            "name": self.name,
            "description": (
                "Run InterProScan domain search against Pfam, TIGRFAM, Gene3D "
                "and other databases via EBI REST API. "
                "Use this when BLAST returns a hypothetical protein, weak e-value "
                "(above 1e-3), or coverage below 50%. "
                "Detects conserved domain signatures even at low sequence identity "
                "where BLAST loses sensitivity. "
                "Do NOT use if BLAST already returned HIGH confidence with a named, "
                "functionally characterized protein. "
                "Takes 30-120 seconds."
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
        Submit InterProScan job, poll for completion, parse top domain hit.

        Never raises: any failure returns AnnotationEvidence with
        ConfidenceTier.UNKNOWN and the error in reasoning.
        """
        sequence = sequence.strip().replace("\n", "").replace(" ", "")

        try:
            job_id = self._submit(sequence)
        except Exception as e:
            return self._failure(f"Failed to submit InterProScan job: {e}")

        try:
            results = self._poll(job_id)
        except Exception as e:
            return self._failure(f"InterProScan polling failed or timed out: {e}")

        try:
            return self._parse(results, sequence)
        except Exception as e:
            return self._failure(f"Failed to parse InterProScan results: {e}")

    def _submit(self, sequence: str) -> str:
        """
        Submit InterProScan job to EBI Job Dispatcher.

        Returns plain text job ID on success.
        """
        response = requests.post(
            _RUN_URL,
            data={
                "email":    self.email,
                "sequence": sequence,
                "goterms":  "false",
                "pathways": "false",
            },
            timeout=30
        )
        response.raise_for_status()
        job_id = response.text.strip()
        if not job_id:
            raise ValueError("EBI returned empty job ID")
        return job_id

    def _poll(self, job_id: str) -> dict:
        """
        Poll EBI status endpoint until FINISHED, then fetch JSON results.

        Status values: RUNNING, FINISHED, ERROR, FAILURE, NOT_FOUND
        """
        elapsed = 0
        while elapsed <= self.max_wait:
            time.sleep(self.poll_interval)
            elapsed += self.poll_interval

            status_resp = requests.get(
                f"{_STATUS_URL}/{job_id}",
                timeout=30
            )
            status_resp.raise_for_status()
            status = status_resp.text.strip()

            if status == "FINISHED":
                result_resp = requests.get(
                    f"{_RESULT_URL}/{job_id}/json",
                    headers={"Accept": "application/json"},
                    timeout=30
                )
                result_resp.raise_for_status()
                return result_resp.json()

            if status in ("ERROR", "FAILURE", "NOT_FOUND"):
                raise RuntimeError(
                    f"InterProScan job {job_id} ended with status: {status}"
                )
            # RUNNING: keep polling

        raise TimeoutError(
            f"InterProScan job {job_id} did not finish within {self.max_wait}s"
        )

    def _parse(self, results: dict, query_sequence: str) -> AnnotationEvidence:
        """
        Parse InterProScan JSON and return the best AnnotationEvidence.

        InterProScan returns matches grouped by sequence. Each match has
        a database source, accession, description, and location entries.
        We select the best hit by database priority then e-value.
        """
        sequences = results.get("results", [])
        if not sequences:
            return self._no_hits()

        # Collect all matches across all databases
        all_matches = []
        for seq_result in sequences:
            for match in seq_result.get("matches", []):
                signature = match.get("signature", {})
                entry = signature.get("entry") or {}
                locations = match.get("locations", [])

                db_name = (
                    signature.get("signatureLibraryRelease", {}).get(
                        "library", ""
                    )
                    or signature.get("library", "")
                    or "unknown"
                )
                accession = signature.get("accession", "unknown")
                description = (
                    entry.get("description")
                    or signature.get("description")
                    or "No description"
                )
                name = (
                    entry.get("name")
                    or signature.get("name")
                    or accession
                )

                # Get best e-value from locations
                evalue = min(
                    (float(loc.get("evalue", 1.0)) for loc in locations),
                    default=1.0
                )

                # Compute coverage from location spans
                coverage = self._compute_coverage(
                    locations, len(query_sequence)
                )

                all_matches.append({
                    "db": db_name,
                    "accession": accession,
                    "name": name,
                    "description": description,
                    "evalue": evalue,
                    "coverage": coverage,
                    "priority": _DB_PRIORITY.get(db_name, 1),
                })

        if not all_matches:
            return self._no_hits()

        # Filter out matches with no meaningful e-value (e-value == 1.0
        # usually means the database doesn't report e-values, e.g. PANTHER)
        scored_matches = [m for m in all_matches if m["evalue"] < 1.0]

        # Fall back to all matches if none have e-values
        candidates = scored_matches if scored_matches else all_matches

        # Select best: lowest e-value first, db priority as tiebreaker
        best = sorted(
            candidates,
            key=lambda m: (m["evalue"], -m["priority"])
        )[0]

        confidence, reasoning = self._assign_confidence(
            best["evalue"], best["coverage"],
            best["description"], best["db"]
        )
        annotation = self._clean_annotation(
            best["description"], best["name"]
        )

        return AnnotationEvidence(
            tool_name=self.name,
            annotation=annotation,
            confidence=confidence,
            score=best["evalue"],
            coverage=best["coverage"],
            hit_id=best["accession"],
            hit_description=best["description"],
            reasoning=reasoning,
            raw_output=str(results)[:500]
        )

    def _compute_coverage(self, locations: list,
                          query_length: int) -> float:
        """Compute fraction of query sequence covered by domain locations."""
        if not locations or query_length <= 0:
            return 0.0

        # Union of all location spans
        covered = set()
        for loc in locations:
            start = int(loc.get("start", 0))
            end   = int(loc.get("end", 0))
            covered.update(range(start, end + 1))

        return len(covered) / query_length

    def _assign_confidence(self, evalue: float, coverage: float,
                           description: str, db: str) -> tuple:
        """
        Map InterProScan hit metrics to ConfidenceTier with reasoning.

        Returns (ConfidenceTier, reasoning_string).
        """
        is_vague = any(
            term in description.lower()
            for term in ["domain of unknown function", "duf",
                         "uncharacterized", "putative", "hypothetical"]
        )

        if evalue <= _HIGH_EVALUE and coverage >= _HIGH_COVERAGE and not is_vague:
            tier = ConfidenceTier.HIGH
            reason = (
                f"Strong domain hit from {db}: e-value {evalue:.2e}, "
                f"{coverage:.0%} query coverage. Domain is well-characterized. "
                f"Annotation supported without further escalation."
            )
        elif evalue <= _MOD_EVALUE and coverage >= _MIN_COVERAGE:
            tier = ConfidenceTier.MODERATE
            if is_vague:
                reason = (
                    f"Domain of unknown function from {db} "
                    f"(e-value {evalue:.2e}, {coverage:.0%} coverage). "
                    f"Signal is real but function is not yet characterized. "
                    f"Structural search recommended."
                )
            else:
                reason = (
                    f"Moderate domain hit from {db}: e-value {evalue:.2e}, "
                    f"{coverage:.0%} coverage. Domain has functional annotation "
                    f"but below high-confidence threshold."
                )
        else:
            tier = ConfidenceTier.LOW
            reason = (
                f"Weak or no domain signal from {db}: e-value {evalue:.2e}, "
                f"{coverage:.0%} coverage. Evidence insufficient. "
                f"Escalate to structural search."
            )

        return tier, reason

    def _clean_annotation(self, description: str, name: str) -> str:
        """Produce a clean annotation label from the InterProScan hit."""
        label = description.strip()
        if not label or label == "No description":
            label = name
        return label[:80] if len(label) > 80 else label

    def _no_hits(self) -> AnnotationEvidence:
        """Return UNKNOWN evidence when InterProScan finds no domain hits."""
        return AnnotationEvidence(
            tool_name=self.name,
            annotation="No domain hits",
            confidence=ConfidenceTier.UNKNOWN,
            reasoning=(
                "InterProScan found no domain hits above threshold across "
                "Pfam, TIGRFAM, Gene3D, and other databases. Sequence has "
                "no detectable conserved domain architecture. Escalate to "
                "structural similarity search via RCSB PDB."
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