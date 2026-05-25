"""
ORACLE command-line entry point.

Wires together the agent, tools, and report renderer into a single
runnable command. Reads a FASTA file, annotates each sequence, and
writes reports to disk.

Usage:
    oracle -i sequences.fasta
    oracle -i sequences.fasta -o reports/ --format json
    oracle -i sequences.fasta --verbose
    oracle -i sequences.fasta --dry-run

Environment variables required:
    ANTHROPIC_API_KEY   — Anthropic API key for the agent loop
    NCBI_API_KEY        — Optional. Raises NCBI rate limit from 3 to 10 req/s

Example .env file (never commit this):
    ANTHROPIC_API_KEY=sk-ant-...
    NCBI_API_KEY=your_ncbi_key_here
"""

import argparse
import os
import sys
from pathlib import Path

from .agent.oracle import OracleAgent
from .tools.registry import ToolRegistry
from .tools.blast import BlastTool
from .tools.hmmer import HmmerTool
from .tools.structure import RCSBStructureTool
from .tools.composition import CompositionTool
from .utils.report import render_text, render_json, save_report


def parse_fasta(fasta_path: str) -> list[tuple[str, str]]:
    """
    Parse a FASTA file into (sequence_id, sequence) tuples.

    Handles multi-sequence FASTA files. Strips whitespace from sequences
    and uses the first word of the header line as the sequence ID.

    Args:
        fasta_path: Path to the FASTA file

    Returns:
        List of (sequence_id, sequence) tuples

    Raises:
        ValueError: If file is empty or contains no valid sequences
        FileNotFoundError: If path does not exist
    """
    path = Path(fasta_path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    sequences = []
    current_id = None
    current_seq = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Save previous sequence if we have one
                if current_id and current_seq:
                    sequences.append((current_id, "".join(current_seq)))
                # Start new sequence, use first word of header as ID
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())

    # Don't forget the last sequence
    if current_id and current_seq:
        sequences.append((current_id, "".join(current_seq)))

    if not sequences:
        raise ValueError(f"No valid sequences found in {fasta_path}")

    return sequences


def build_registry() -> ToolRegistry:
    """
    Build and populate the tool registry with all available tools.

    Add new tools here as they are implemented. The agent will
    automatically receive their schemas and can choose to call them.
    """
    registry = ToolRegistry()
    registry.register(BlastTool())
    registry.register(HmmerTool())
    registry.register(RCSBStructureTool())
    registry.register(CompositionTool())

    return registry


def check_environment() -> list[str]:
    """
    Check required environment variables are set.

    Returns list of missing variables. Empty list means all good.
    """
    missing = []
    if not os.environ.get("ANTHROPIC_API_KEY"):
        missing.append("ANTHROPIC_API_KEY")
    return missing


def main():
    parser = argparse.ArgumentParser(
        prog="oracle",
        description=(
            "ORACLE: Autonomous annotation agent for uncharacterized protein sequences.\n"
            "Uses confidence-layered evidence reasoning across homology, domain,\n"
            "structural, and literature sources."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        metavar="FASTA",
        help="Input FASTA file. Single or multi-sequence."
    )
    parser.add_argument(
        "-o", "--output",
        default="reports",
        metavar="DIR",
        help="Output directory for reports. Default: reports/"
    )
    parser.add_argument(
        "--format",
        choices=["text", "json", "both"],
        default="both",
        help="Output format. Default: both"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print each agent step to stdout as it runs."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Parse input and validate environment without running the agent."
    )

    args = parser.parse_args()

    # Load .env file if present (simple implementation, no dotenv dependency)
    env_file = Path(".env")
    if env_file.exists():
        with open(env_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#") and "=" in line:
                    key, _, value = line.partition("=")
                    os.environ.setdefault(key.strip(), value.strip())

    # Check environment before doing anything else
    missing = check_environment()
    if missing:
        print(f"Error: missing required environment variables: {', '.join(missing)}")
        print("Set them in your shell or in a .env file in the project root.")
        sys.exit(1)

    # Parse input file
    try:
        sequences = parse_fasta(args.input)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error reading input: {e}")
        sys.exit(1)

    print(f"ORACLE: found {len(sequences)} sequence(s) in {args.input}")

    if args.dry_run:
        print("Dry run complete. Environment and input look good.")
        for seq_id, seq in sequences:
            print(f"  {seq_id}: {len(seq)} aa")
        sys.exit(0)

    # Build registry and agent
    registry = build_registry()
    agent = OracleAgent(registry, verbose=args.verbose)

    print(f"Tools registered: {registry.available_tools()}")
    print(f"Output directory: {args.output}")
    print()

    # Annotate each sequence
    for i, (seq_id, seq) in enumerate(sequences, 1):
        print(f"[{i}/{len(sequences)}] Annotating {seq_id} ({len(seq)} aa)...")

        try:
            report = agent.annotate(seq_id, seq)
        except Exception as e:
            print(f"  Error annotating {seq_id}: {e}")
            continue

        # Print text report to terminal
        if args.verbose or len(sequences) == 1:
            print(render_text(report))

        # Save to disk
        written = save_report(report, output_dir=args.output, fmt=args.format)
        for path in written:
            print(f"  Saved: {path}")

        print(f"  Result: {report.final_annotation} [{report.final_confidence.name}]")
        print()

    print("ORACLE: annotation complete.")


if __name__ == "__main__":
    main()