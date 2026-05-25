"""
ORACLE agent loop.

This is the core of the system. The agent receives a protein sequence,
reasons about which tools to call using the Anthropic tool use API,
dispatches each call through the ToolRegistry, and builds an OracleReport
from the accumulated evidence.

The loop follows the Anthropic tool use pattern:
  1. Send message with tool schemas to the API
  2. Model responds with tool_use blocks (which tool to call and with what)
  3. Execute the tool via registry, get AnnotationEvidence back
  4. Send tool result back to the model
  5. Repeat until model responds with text only (done reasoning)

The agent never calls tools directly. All dispatch goes through the
registry so tools are swappable without changing agent logic.
"""

import json
import os
from dataclasses import asdict

import anthropic

from .prompts import SYSTEM_PROMPT, build_user_message
from ..models import (
    AnnotationEvidence, AgentStep, OracleReport, ConfidenceTier
)
from ..tools.registry import ToolRegistry


# Maximum tool calls per run. Prevents runaway loops on unexpected model
# behavior. With 3 tools in the ladder this should never be hit in normal
# operation. Raise if you add more tools.
_MAX_TOOL_CALLS = 6

# Model to use for the agent loop
_MODEL = "claude-sonnet-4-5-20250929"


class OracleAgent:
    """
    Autonomous protein annotation agent.

    Wraps the Anthropic tool use loop and manages the evidence chain.
    Each call to annotate() is a complete independent annotation run.

    Usage:
        from oracle.tools.blast import BlastTool
        from oracle.tools.hmmer import HmmerTool

        registry = ToolRegistry()
        registry.register(BlastTool())
        registry.register(HmmerTool())

        agent = OracleAgent(registry)
        report = agent.annotate("seq_001", "MKTIIALSYIFCLVFA...")
        print(report.final_annotation)
    """

    def __init__(self, registry: ToolRegistry, verbose: bool = False):
        """
        Args:
            registry: Populated ToolRegistry with all available tools.
            verbose:  If True, print each agent step to stdout as it happens.
                      Useful during development; turn off for production.
        """
        self.registry = registry
        self.verbose = verbose
        self.client = anthropic.Anthropic(
            api_key=os.environ.get("ANTHROPIC_API_KEY")
        )

    def annotate(self, sequence_id: str, sequence: str) -> OracleReport:
        """
        Run the full annotation loop for one protein sequence.

        Args:
            sequence_id: Identifier for this sequence (FASTA header, accession)
            sequence:    Clean amino acid sequence string, no header or spaces

        Returns:
            OracleReport with final annotation, confidence tier, and full
            evidence chain showing every tool call and reasoning step.
        """
        sequence = sequence.strip().replace("\n", "").replace(" ", "")

        # Track state across the tool use loop
        messages = [
            {"role": "user", "content": build_user_message(sequence_id, sequence)}
        ]
        evidence_chain: list[AgentStep] = []
        skipped_tools: list[dict] = []
        tool_call_count = 0
        tools_called: set[str] = set()

        if self.verbose:
            print(f"\n[ORACLE] Annotating {sequence_id} ({len(sequence)} aa)")

        # Main agent loop
        while tool_call_count < _MAX_TOOL_CALLS:

            response = self.client.messages.create(
                model=_MODEL,
                max_tokens=2048,
                system=SYSTEM_PROMPT,
                tools=self.registry.get_schemas(),
                messages=messages,
            )

            # Append assistant response to message history
            messages.append({"role": "assistant", "content": response.content})

            # Check stop condition: model responded with text only, done reasoning
            if response.stop_reason == "end_turn":
                if self.verbose:
                    print("[ORACLE] Agent finished reasoning")
                break

            # Process all tool use blocks in this response
            # (model may request multiple tools in one turn, though unlikely
            # given our system prompt instructs one tool at a time)
            tool_results = []
            for block in response.content:
                if block.type != "tool_use":
                    continue

                tool_name = block.name
                tool_input = block.input
                tool_call_count += 1

                # Guard against the model calling the same tool twice
                if tool_name in tools_called:
                    if self.verbose:
                        print(f"[ORACLE] Skipping duplicate call to {tool_name}")
                    skipped_tools.append({
                        "tool": tool_name,
                        "reason": "Agent attempted duplicate call, skipped."
                    })
                    # Return a fake result so the conversation stays valid
                    tool_results.append(self._format_tool_result(
                        block.id, tool_name,
                        AnnotationEvidence(
                            tool_name=tool_name,
                            annotation="Duplicate call skipped",
                            confidence=ConfidenceTier.UNKNOWN,
                            reasoning="This tool was already called this run."
                        )
                    ))
                    continue

                tools_called.add(tool_name)

                if self.verbose:
                    print(f"[ORACLE] Calling tool: {tool_name}")

                # Execute tool via registry
                evidence = self.registry.run(
                    tool_name,
                    sequence=tool_input.get("sequence", sequence),
                    **{k: v for k, v in tool_input.items() if k != "sequence"}
                )

                if self.verbose:
                    print(f"[ORACLE] {tool_name} -> {evidence.confidence.name}: "
                          f"{evidence.annotation}")
                    if evidence.annotation in ("Tool failure", "Tool not found"):
                        print(f"[ORACLE] Tool reasoning: {evidence.reasoning}")

                # Extract the agent's stated reason for this tool call
                # from any preceding text block in the response
                reason_for_choice = self._extract_reasoning(response.content)

                # Record this step in the evidence chain.
                # escalate=False when confidence is HIGH (no need to continue)
                # or when this is the last registered tool (nothing left to call)
                available = self.registry.available_tools()
                is_last_tool = (tool_name == available[-1])
                should_escalate = (
                    evidence.confidence.value < ConfidenceTier.HIGH.value
                    and not is_last_tool
                )
                step = AgentStep(
                    step_number=len(evidence_chain) + 1,
                    tool_chosen=tool_name,
                    reason_for_choice=reason_for_choice,
                    evidence=evidence,
                    updated_hypothesis=evidence.annotation,
                    escalate=should_escalate
                )
                evidence_chain.append(step)

                tool_results.append(
                    self._format_tool_result(block.id, tool_name, evidence)
                )

            # Send all tool results back to the model in one message
            if tool_results:
                messages.append({
                    "role": "user",
                    "content": tool_results
                })

        # Extract final annotation from the last text response
        final_text = self._extract_final_text(messages)
        final_annotation, final_confidence, warnings = self._parse_final_response(
            final_text, evidence_chain
        )

        # Note any registered tools that were never called
        for tool_name in self.registry.available_tools():
            if tool_name not in tools_called:
                skipped_tools.append({
                    "tool": tool_name,
                    "reason": "Not needed given evidence from earlier tools."
                })

        return OracleReport(
            sequence_id=sequence_id,
            sequence_length=len(sequence),
            final_annotation=final_annotation,
            final_confidence=final_confidence,
            evidence_chain=evidence_chain,
            skipped_tools=skipped_tools,
            summary_reasoning=final_text,
            warnings=warnings,
        )

    def _format_tool_result(
        self, tool_use_id: str, tool_name: str, evidence: AnnotationEvidence
    ) -> dict:
        """
        Format an AnnotationEvidence object as an Anthropic tool result block.

        The model receives this as structured JSON so it can reason about
        specific fields like confidence, score, and coverage rather than
        parsing free text.
        """
        # Convert dataclass to dict, handling the ConfidenceTier enum
        evidence_dict = {
            "tool_name": evidence.tool_name,
            "annotation": evidence.annotation,
            "confidence": evidence.confidence.name,
            "score": evidence.score,
            "coverage": evidence.coverage,
            "hit_id": evidence.hit_id,
            "hit_description": evidence.hit_description,
            "reasoning": evidence.reasoning,
        }

        return {
            "type": "tool_result",
            "tool_use_id": tool_use_id,
            "content": json.dumps(evidence_dict, indent=2)
        }

    def _extract_reasoning(self, content_blocks: list) -> str:
        """
        Pull any text the model wrote before its tool call.

        The agent often explains its reasoning in a text block before
        requesting a tool. We capture this as the reason_for_choice
        in AgentStep so the evidence chain is readable.
        """
        for block in content_blocks:
            if hasattr(block, "type") and block.type == "text" and block.text.strip():
                return block.text.strip()
        return "No explicit reasoning provided before tool call."

    def _extract_final_text(self, messages: list) -> str:
        """
        Find the last assistant text response in the message history.

        This is the agent's final summary after all tool calls are complete.
        """
        for message in reversed(messages):
            if message.get("role") != "assistant":
                continue
            for block in message.get("content", []):
                if hasattr(block, "type") and block.type == "text":
                    return block.text.strip()
                # Handle dict format as well as object format
                if isinstance(block, dict) and block.get("type") == "text":
                    return block.get("text", "").strip()
        return "No final summary produced."

    def _parse_final_response(
        self, final_text: str, evidence_chain: list[AgentStep]
    ) -> tuple[str, ConfidenceTier, list[str]]:
        """
        Extract structured fields from the agent's final text response.

        Derives final annotation and confidence from the evidence chain.
        Uses the highest confidence result across all steps rather than
        just the last step, since the agent may stop early on a strong
        hit before running all tools.

        Returns: (final_annotation, final_confidence, warnings)
        """
        warnings = []

        if not evidence_chain:
            return "Uncharacterized protein", ConfidenceTier.UNKNOWN, [
                "No tools were called. Sequence may be invalid or agent loop failed."
            ]

        # Find the best confidence result across all steps, excluding
        # tool failures and duplicate call artifacts which are not real evidence.
        _failure_annotations = {
            "tool failure", "tool not found", "duplicate call skipped"
        }
        valid_steps = [
            s for s in evidence_chain
            if s.evidence.annotation.lower() not in _failure_annotations
        ]

        # If all steps failed fall back to the full chain
        if not valid_steps:
            valid_steps = evidence_chain

        best_step = max(
            valid_steps,
            key=lambda s: s.evidence.confidence.value
        )
        final_confidence = best_step.evidence.confidence

        # For HIGH confidence results, use the best step's annotation.
        # For anything below HIGH, use the LAST valid step's annotation
        # since it reflects the most complete evidence. Then enrich it
        # with BLAST context if available, to distinguish:
        #   A) "Uncharacterized phage protein (homolog in Vicingus sp.)"
        #      — BLAST found something, just uninformative
        #   B) "Uncharacterized protein (no database representatives)"
        #      — genuinely orphan, nothing anywhere
        if final_confidence == ConfidenceTier.HIGH:
            final_annotation = best_step.evidence.annotation
        else:
            last_step = valid_steps[-1]
            final_annotation = last_step.evidence.annotation
            final_confidence = last_step.evidence.confidence

            # Check if BLAST found a hit even if uninformative
            blast_step = next(
                (s for s in valid_steps if s.tool_chosen == "blast_nr"),
                None
            )
            blast_had_hit = (
                blast_step is not None
                and blast_step.evidence.hit_id is not None
                and blast_step.evidence.confidence != ConfidenceTier.UNKNOWN
                and "no significant" not in blast_step.evidence.annotation.lower()
                and "tool failure" not in blast_step.evidence.annotation.lower()
            )

            # If the final annotation is a composition or no-hits result,
            # add BLAST context to differentiate the two dark matter scenarios
            uninformative_finals = {
                "no pdb hits", "no domain hits", "no structural hits",
                "tool not available"
            }
            final_lower = final_annotation.lower()
            is_composition_result = last_step.tool_chosen == "composition"
            is_no_hits = any(t in final_lower for t in uninformative_finals)

            if is_composition_result or is_no_hits:
                if blast_had_hit and blast_step.evidence.hit_id:
                    # Scenario A: conserved but uncharacterized
                    organism = ""
                    hit_desc = blast_step.evidence.hit_description or ""
                    if "[" in hit_desc:
                        organism = hit_desc[hit_desc.index("[")+1:].split("]")[0]
                    if organism:
                        final_annotation = (
                            f"Uncharacterized protein "
                            f"(homolog detected in {organism})"
                        )
                    else:
                        final_annotation = (
                            "Uncharacterized protein "
                            "(conserved homolog in database, function unknown)"
                        )
                else:
                    # Scenario B: true orphan
                    final_annotation = (
                        "Uncharacterized protein "
                        "(no database representatives found)"
                    )

        # If we still have a failure annotation, the agent's summary
        # reasoning is more informative than "Tool failure". Extract
        # the putative annotation from the summary text if present.
        if final_annotation.lower() in _failure_annotations and final_text:
            # Look for "Final Annotation:" pattern in agent summary
            for line in final_text.splitlines():
                line = line.strip().lstrip("*# ")
                if line.lower().startswith("final annotation"):
                    parts = line.split(":", 1)
                    if len(parts) > 1:
                        candidate = parts[1].strip().lstrip("*# ")
                        if candidate:
                            final_annotation = candidate
                            final_confidence = ConfidenceTier.LOW
                            break

        # If the agent's summary text explicitly mentions HIGH confidence
        # and we only have MODERATE from tools, trust the agent's reasoning.
        # This handles cases where the hit description is parsed as hypothetical
        # but the full evidence (coverage, e-value, context) is clearly strong.
        if final_confidence == ConfidenceTier.MODERATE:
            summary_upper = final_text.upper()
            if "HIGH" in summary_upper and "CONFIDENCE" in summary_upper:
                final_confidence = ConfidenceTier.HIGH

        # Check for conflicts between tools - only consider steps that
        # produced real functional annotations, not absence-of-evidence results
        _non_annotations = {
            "hypothetical", "no significant", "no domain", "no pdb",
            "tool failure", "consensus disorder", "uncharacterized",
            "tool not available", "no structural", "no hits",
            "coil", "signal peptide", "transmembrane"
        }
        annotations = [
            step.evidence.annotation for step in evidence_chain
            if not any(t in step.evidence.annotation.lower() for t in _non_annotations)
            and step.evidence.annotation.lower() not in _failure_annotations
        ]
        unique_annotations = set(annotations)
        if len(unique_annotations) > 1:
            warnings.append(
                f"Conflicting annotations across tools: {unique_annotations}. "
                f"Experimental validation recommended."
            )

        # Flag if we ended on a weak result
        if final_confidence in (ConfidenceTier.LOW, ConfidenceTier.UNKNOWN):
            warnings.append(
                "Annotation is low confidence. Sequence may be genuine dark matter "
                "with no current database representatives."
            )

        return final_annotation, final_confidence, warnings