"""
Base interface for all ORACLE annotation tools.

Every tool inherits from AnnotationTool and implements run().
This consistent interface lets the agent treat all tools the same
way regardless of whether they call a local binary, a REST API,
or the Anthropic API.

The tool_schema property is what gets passed to the Anthropic API
as a tool definition. The description field in the schema is the
most important part: it tells the model when to use this tool
and when not to. Vague descriptions produce poor agent decisions.
"""

from abc import ABC, abstractmethod
from ..models import AnnotationEvidence


class AnnotationTool(ABC):
    """
    Abstract base class for all ORACLE tools.

    Subclass this for every tool the agent can use. The agent
    never calls tools directly; it selects them by name and the
    tool registry handles dispatch.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """
        Unique tool identifier used in agent tool calls and report output.
        Keep this lowercase with underscores: 'blast_nr', 'hmmer_pfam'.
        """
        ...

    @property
    @abstractmethod
    def tool_schema(self) -> dict:
        """
        Anthropic API tool definition schema.

        Structure:
        {
            "name": self.name,
            "description": "When to use this tool and what it does.
                           Be specific about when NOT to use it too.",
            "input_schema": {
                "type": "object",
                "properties": { ... },
                "required": [ ... ]
            }
        }
        """
        ...

    @abstractmethod
    def run(self, sequence: str, **kwargs) -> AnnotationEvidence:
        """
        Execute the tool on the given sequence.

        Args:
            sequence: Amino acid or nucleotide sequence string (no FASTA header)
            **kwargs: Tool-specific parameters from the agent's tool call

        Returns:
            AnnotationEvidence with the tool's best interpretation.
            Never raises on tool failure; return AnnotationEvidence with
            ConfidenceTier.UNKNOWN and explain the failure in reasoning.
        """
        ...

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(name='{self.name}')"