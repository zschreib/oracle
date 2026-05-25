"""
Tool registry for ORACLE.

The registry is the single source of truth for which tools are available.
The agent asks the registry for tool schemas (to pass to Anthropic API)
and dispatches tool calls through it. Adding a new tool means registering
it here; nothing else needs to change.
"""

from .base import AnnotationTool
from ..models import AnnotationEvidence, ConfidenceTier


class ToolRegistry:
    """
    Manages available annotation tools and handles agent dispatch.

    Usage:
        registry = ToolRegistry()
        registry.register(BlastTool())
        registry.register(HmmerTool())

        # Get schemas for Anthropic API
        schemas = registry.get_schemas()

        # Dispatch a tool call from the agent
        evidence = registry.run("blast_nr", sequence="MKTII...")
    """

    def __init__(self):
        # tool name -> tool instance
        self._tools: dict[str, AnnotationTool] = {}

    def register(self, tool: AnnotationTool) -> None:
        """Register a tool. Overwrites if name already exists."""
        self._tools[tool.name] = tool

    def get_schemas(self) -> list[dict]:
        """
        Return all tool schemas for the Anthropic API tools parameter.
        Called once at agent initialization.
        """
        return [tool.tool_schema for tool in self._tools.values()]

    def run(self, tool_name: str, sequence: str, **kwargs) -> AnnotationEvidence:
        """
        Dispatch a tool call from the agent.

        Args:
            tool_name: Must match a registered tool's name property
            sequence:  Query sequence string
            **kwargs:  Additional parameters from the agent's tool call

        Returns:
            AnnotationEvidence from the tool, or a UNKNOWN-tier evidence
            object if the tool name is not found (agent gets honest feedback
            rather than a crash).
        """
        if tool_name not in self._tools:
            # Return honest failure rather than raising so the agent
            # can reason about the missing tool and proceed
            return AnnotationEvidence(
                tool_name=tool_name,
                annotation="Tool not found",
                confidence=ConfidenceTier.UNKNOWN,
                reasoning=f"Tool '{tool_name}' is not registered. "
                          f"Available tools: {list(self._tools.keys())}"
            )

        return self._tools[tool_name].run(sequence, **kwargs)

    def available_tools(self) -> list[str]:
        """Return list of registered tool names."""
        return list(self._tools.keys())

    def __len__(self) -> int:
        return len(self._tools)