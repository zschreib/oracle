"""
Report renderer for ORACLE.

Converts an OracleReport into human-readable text and machine-readable
JSON. The text format is designed to be readable in a terminal and
copyable into a lab notebook. The JSON format preserves all fields
for downstream processing or database ingestion.

Both formats include the full evidence chain so the annotation is
always traceable back to the specific tool outputs that produced it.
"""

import json
from datetime import datetime
from ..models import OracleReport, ConfidenceTier


# Visual markers for confidence tiers in terminal output
_TIER_MARKERS = {
    ConfidenceTier.HIGH:     "████ HIGH",
    ConfidenceTier.MODERATE: "███░ MODERATE",
    ConfidenceTier.LOW:      "██░░ LOW",
    ConfidenceTier.UNKNOWN:  "█░░░ UNKNOWN",
}


def render_text(report: OracleReport) -> str:
    """
    Render an OracleReport as a formatted text report.

    Designed to be readable in a terminal and meaningful without
    any additional context. Every annotation decision is traceable
    through the evidence chain section.
    """
    lines = []
    width = 68
    bar = "=" * width

    # Header
    lines.append(bar)
    lines.append("ORACLE ANNOTATION REPORT")
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(bar)

    # Sequence info
    lines.append("")
    lines.append("SEQUENCE")
    lines.append(f"  ID:     {report.sequence_id}")
    lines.append(f"  Length: {report.sequence_length} amino acids")

    # Final annotation (most important, shown prominently)
    lines.append("")
    lines.append("ANNOTATION")
    lines.append(f"  {report.final_annotation}")
    lines.append("")
    lines.append("CONFIDENCE")
    lines.append(f"  {_TIER_MARKERS[report.final_confidence]}")

    # Warnings
    if report.warnings:
        lines.append("")
        lines.append("WARNINGS")
        for warning in report.warnings:
            lines.append(f"  ! {warning}")

    # Evidence chain - the core of what makes ORACLE different from a pipeline
    lines.append("")
    lines.append("EVIDENCE CHAIN")
    lines.append("-" * width)

    if not report.evidence_chain:
        lines.append("  No tool calls recorded.")
    else:
        for step in report.evidence_chain:
            ev = step.evidence
            lines.append(f"  Step {step.step_number}: {step.tool_chosen.upper()}")
            lines.append(f"  Reason called: {step.reason_for_choice}")
            lines.append(f"  Result:        {ev.annotation}")
            lines.append(f"  Confidence:    {ev.confidence.name}")

            # Show quantitative metrics when available
            if ev.score is not None:
                lines.append(f"  Score:         {ev.score:.2e}")
            if ev.coverage is not None:
                lines.append(f"  Coverage:      {ev.coverage:.0%}")
            if ev.hit_id:
                lines.append(f"  Hit ID:        {ev.hit_id}")

            lines.append(f"  Reasoning:     {ev.reasoning}")
            escalate_str = "Yes, escalate to next tool" if step.escalate else "No, evidence sufficient"
            lines.append(f"  Escalate:      {escalate_str}")
            lines.append("")

    # Skipped tools
    if report.skipped_tools:
        lines.append("TOOLS NOT CALLED")
        lines.append("-" * width)
        for skip in report.skipped_tools:
            lines.append(f"  {skip['tool']}: {skip['reason']}")
        lines.append("")

    # Agent summary reasoning
    if report.summary_reasoning and report.summary_reasoning != "No final summary produced.":
        lines.append("AGENT SUMMARY")
        lines.append("-" * width)
        # Wrap long lines for terminal readability
        for line in report.summary_reasoning.splitlines():
            if len(line) <= width:
                lines.append(f"  {line}")
            else:
                # Simple word wrap at width characters
                words = line.split()
                current = "  "
                for word in words:
                    if len(current) + len(word) + 1 > width:
                        lines.append(current)
                        current = f"  {word}"
                    else:
                        current = f"{current} {word}" if current.strip() else f"  {word}"
                if current.strip():
                    lines.append(current)
        lines.append("")

    lines.append(bar)
    return "\n".join(lines)


def render_json(report: OracleReport) -> str:
    """
    Render an OracleReport as JSON for downstream processing.

    All fields are included. ConfidenceTier enums are serialized as
    their string names (HIGH, MODERATE, LOW, UNKNOWN) for readability.
    """
    def _serialize(obj):
        """Handle types that json.dumps can't serialize natively."""
        if isinstance(obj, ConfidenceTier):
            return obj.name
        raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

    report_dict = {
        "oracle_version": "1.0.0",
        "generated": datetime.now().isoformat(),
        "sequence_id": report.sequence_id,
        "sequence_length": report.sequence_length,
        "final_annotation": report.final_annotation,
        "final_confidence": report.final_confidence.name,
        "warnings": report.warnings,
        "evidence_chain": [
            {
                "step": step.step_number,
                "tool": step.tool_chosen,
                "reason_for_choice": step.reason_for_choice,
                "annotation": step.evidence.annotation,
                "confidence": step.evidence.confidence.name,
                "score": step.evidence.score,
                "coverage": step.evidence.coverage,
                "hit_id": step.evidence.hit_id,
                "hit_description": step.evidence.hit_description,
                "reasoning": step.evidence.reasoning,
                "escalate": step.escalate,
            }
            for step in report.evidence_chain
        ],
        "skipped_tools": report.skipped_tools,
        "summary_reasoning": report.summary_reasoning,
    }

    return json.dumps(report_dict, default=_serialize, indent=2)


def save_report(report: OracleReport, output_dir: str = "reports",
                fmt: str = "both") -> list[str]:
    """
    Save report to disk in text, JSON, or both formats.

    Args:
        report:     OracleReport to save
        output_dir: Directory to write files into (created if absent)
        fmt:        'text', 'json', or 'both'

    Returns:
        List of file paths written.
    """
    import os
    os.makedirs(output_dir, exist_ok=True)

    # Sanitize sequence ID for use as filename.
    # Windows disallows: \ / : * ? " < > | and spaces.
    import re
    safe_id = re.sub(r'[\\/:*?"<>|\s]', "_", report.sequence_id)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base = os.path.join(output_dir, f"{safe_id}_{timestamp}")

    written = []

    if fmt in ("text", "both"):
        path = f"{base}.report.txt"
        with open(path, "w", encoding="utf-8") as f:
            f.write(render_text(report))
        written.append(path)

    if fmt in ("json", "both"):
        path = f"{base}.report.json"
        with open(path, "w", encoding="utf-8") as f:
            f.write(render_json(report))
        written.append(path)

    return written