"""
Core data models for ORACLE.

These dataclasses define the shared vocabulary between the agent,
tools, and report generator. Every tool returns an AnnotationEvidence
object so the agent reasons over a consistent interface regardless of
which tool produced the result.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


class ConfidenceTier(Enum):
    """
    Tiered confidence system for annotations.

    Tiers reflect the strength and type of evidence, not just a score.
    The agent uses these to decide whether to escalate to additional tools.

    HIGH:     Strong homology hit (e-value < 1e-10, coverage > 70%) or
              multiple independent lines of evidence agree.
    MODERATE: Domain hit without full sequence homology, or structural
              match with functional support from literature.
    LOW:      Weak or partial evidence from a single source. Annotation
              is a hypothesis, not an assignment.
    UNKNOWN:  All tools exhausted with no interpretable signal. Sequence
              is genuinely dark matter by current databases.
    """
    HIGH = 4
    MODERATE = 3
    LOW = 2
    UNKNOWN = 1


@dataclass
class AnnotationEvidence:
    """
    Structured output from a single tool run.

    Each tool populates this with its best interpretation so the agent
    can reason over results without parsing raw tool output directly.
    The reasoning field is critical: it captures why the tool produced
    this result, not just what it found.
    """
    tool_name: str
    annotation: str                        # human-readable functional label
    confidence: ConfidenceTier
    score: Optional[float] = None          # e-value, bitscore, TM-score etc
    coverage: Optional[float] = None       # fraction of query covered by hit
    hit_id: Optional[str] = None          # accession or identifier of top hit
    hit_description: Optional[str] = None # full description of top hit
    reasoning: str = ""                    # why this annotation was assigned
    raw_output: Optional[str] = None      # raw tool output for debugging


@dataclass
class AgentStep:
    """
    Records a single decision step in the agent reasoning loop.

    The agent appends one of these for every tool it decides to run,
    including its reasoning for running that tool and its interpretation
    of the result. Together these form the evidence chain in the report.
    """
    step_number: int
    tool_chosen: str
    reason_for_choice: str       # why the agent chose this tool at this step
    evidence: AnnotationEvidence
    updated_hypothesis: str      # agent's working annotation after this step
    escalate: bool               # whether agent decided to run another tool


@dataclass
class OracleReport:
    """
    Final output of a complete ORACLE annotation run.

    Contains the final annotation, confidence tier, full evidence chain,
    and a list of tools that were skipped and why. The skipped_tools
    field is important: it shows the agent made deliberate decisions
    about tool use rather than running everything blindly.
    """
    sequence_id: str
    sequence_length: int
    final_annotation: str
    final_confidence: ConfidenceTier
    evidence_chain: list[AgentStep] = field(default_factory=list)
    skipped_tools: list[dict] = field(default_factory=list)  # {tool, reason}
    summary_reasoning: str = ""   # agent's narrative explanation of the result
    warnings: list[str] = field(default_factory=list)  # conflicts or caveats