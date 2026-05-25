"""
Mock tool for ORACLE development and testing.

Simulates the full tool escalation ladder without making any API calls.
The mock returns pre-configured AnnotationEvidence objects based on
the sequence input, letting you test every agent decision path cheaply.

Usage in tests:
    from oracle.tools.mock import MockTool
    from oracle.models import ConfidenceTier

    # Simulate a strong BLAST hit - agent should stop here
    tool = MockTool(
        name="blast_nr",
        confidence=ConfidenceTier.HIGH,
        annotation="tail fiber protein",
        score=1e-50,
        coverage=0.92
    )

Usage in development:
    # Simulate dark matter sequence - agent should escalate all the way
    blast_mock = MockTool("blast_nr", ConfidenceTier.LOW, "No significant hits")
    hmmer_mock = MockTool("hmmer_pfam", ConfidenceTier.LOW, "No domain hits")
    foldseek_mock = MockTool("rcsb_pdb", ConfidenceTier.MODERATE, "Putative hydrolase fold")
"""

from oracle.tools.base import AnnotationTool
from oracle.models import AnnotationEvidence, ConfidenceTier


class MockTool(AnnotationTool):
    """
    Configurable mock that returns fixed AnnotationEvidence.

    Useful for testing agent escalation logic without API calls.
    Set delay_seconds > 0 to simulate realistic tool latency in
    integration tests.
    """

    def __init__(
        self,
        name: str,
        confidence: ConfidenceTier,
        annotation: str,
        score: float = None,
        coverage: float = None,
        hit_id: str = "mock_hit_001",
        hit_description: str = None,
        reasoning: str = None,
        delay_seconds: float = 0.0,
    ):
        self._name = name
        self._confidence = confidence
        self._annotation = annotation
        self._score = score
        self._coverage = coverage
        self._hit_id = hit_id
        self._hit_description = hit_description or annotation
        self._reasoning = reasoning or self._default_reasoning()
        self._delay_seconds = delay_seconds

        # Track calls for test assertions
        self.call_count = 0
        self.last_sequence = None

    @property
    def name(self) -> str:
        return self._name

    @property
    def tool_schema(self) -> dict:
        return {
            "name": self._name,
            "description": (
                f"Mock tool simulating {self._name}. "
                f"Always returns {self._confidence.name} confidence. "
                f"For testing only."
            ),
            "input_schema": {
                "type": "object",
                "properties": {
                    "sequence": {
                        "type": "string",
                        "description": "Amino acid sequence string."
                    }
                },
                "required": ["sequence"]
            }
        }

    def run(self, sequence: str, **kwargs) -> AnnotationEvidence:
        """
        Return the pre-configured evidence immediately.

        Records call count and sequence for test assertions.
        """
        import time

        self.call_count += 1
        self.last_sequence = sequence

        if self._delay_seconds > 0:
            time.sleep(self._delay_seconds)

        return AnnotationEvidence(
            tool_name=self._name,
            annotation=self._annotation,
            confidence=self._confidence,
            score=self._score,
            coverage=self._coverage,
            hit_id=self._hit_id,
            hit_description=self._hit_description,
            reasoning=self._reasoning,
            raw_output="[mock output]"
        )

    def _default_reasoning(self) -> str:
        """Generate sensible default reasoning text for each confidence tier."""
        defaults = {
            ConfidenceTier.HIGH: (
                f"Strong homology detected. Mock e-value 1e-50, 92% coverage. "
                f"Annotation '{self._annotation}' is well-supported."
            ),
            ConfidenceTier.MODERATE: (
                f"Moderate evidence for '{self._annotation}'. "
                f"Coverage or e-value below high-confidence threshold. "
                f"Consider escalating to next tool."
            ),
            ConfidenceTier.LOW: (
                f"Weak signal. Annotation '{self._annotation}' is speculative. "
                f"Escalate to next tool in the evidence ladder."
            ),
            ConfidenceTier.UNKNOWN: (
                f"No interpretable signal from {self._name}. "
                f"Tool exhausted without meaningful result."
            ),
        }
        return defaults.get(self._confidence, "Mock reasoning.")


# Pre-built scenario sets for common test cases
# Import these directly in tests rather than constructing mocks manually

def scenario_strong_blast_hit() -> list:
    """Agent should annotate on BLAST alone and skip remaining tools."""
    return [
        MockTool("blast_nr", ConfidenceTier.HIGH,
                 "tail fiber protein", score=1e-50, coverage=0.92),
        MockTool("hmmer_pfam", ConfidenceTier.HIGH,
                 "tail fiber protein", score=250.0, coverage=0.88),
        MockTool("rcsb_pdb", ConfidenceTier.HIGH,
                 "tail fiber protein"),
    ]


def scenario_blast_fails_hmmer_succeeds() -> list:
    """Agent should escalate from BLAST to HMMer and stop there."""
    return [
        MockTool("blast_nr", ConfidenceTier.LOW,
                 "No significant hits", score=0.8, coverage=0.15),
        MockTool("hmmer_pfam", ConfidenceTier.MODERATE,
                 "Lysozyme domain", score=145.0, coverage=0.75),
        MockTool("rcsb_pdb", ConfidenceTier.MODERATE,
                 "Lysozyme fold"),
    ]


def scenario_dark_matter() -> list:
    """All tools fail. Agent should return UNKNOWN with full evidence trail."""
    return [
        MockTool("blast_nr", ConfidenceTier.LOW,
                 "No significant hits", score=1.2, coverage=0.08),
        MockTool("hmmer_pfam", ConfidenceTier.UNKNOWN,
                 "No domain hits"),
        MockTool("rcsb_pdb", ConfidenceTier.LOW,
                 "Weak structural similarity", score=0.42, coverage=0.35),
    ]


def scenario_conflicting_evidence() -> list:
    """BLAST and HMMer disagree. Agent should flag conflict in report."""
    return [
        MockTool("blast_nr", ConfidenceTier.MODERATE,
                 "hypothetical protein",
                 reasoning="Hit is hypothetical. Structural similarity to "
                           "lysozyme family but no named function."),
        MockTool("hmmer_pfam", ConfidenceTier.MODERATE,
                 "Tail spike protein",
                 reasoning="InterProScan domain hit suggests tail spike, "
                           "inconsistent with BLAST lysozyme signal."),
        MockTool("rcsb_pdb", ConfidenceTier.MODERATE,
                 "Lysozyme-like fold"),
    ]