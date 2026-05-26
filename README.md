# ORACLE
**ORf Annotation via Reasoning and Confidence-Layered Evidence**

An autonomous protein annotation agent that combines Anthropic's tool use API with a four-tier evidence ladder to annotate uncharacterized or uninformative protein sequences. ORACLE is built for metagenomic and viral dark matter: proteins with little to no sequence homologs, conserved domains, and/or structural relatives in any database.

---
## Why ORACLE

Automated annotation pipelines run a fixed set of tools and report whatever they find. ORACLE reasons about what it finds. It decides when to escalate, when to stop, and how to interpret conflicting evidence across tools. For well-characterized proteins it stops at BLAST. For genuine dark matter it runs all four tiers and produces a structured hypothesis from sequence composition.

The output is not just an annotation, it is an evidence chain with explicit confidence, reasoning at each step, and specific recommendations for experimental follow-up.

## Evidence ladder

ORACLE runs tools in order, stopping as soon as it has high-confidence evidence. Each tier is only called when the previous one fails to produce a confident functional annotation.

| Tier | Tool | Database | Fires when |
|------|------|----------|------------|
| 1 | BLASTp | NCBI nr | Always first |
| 2 | InterProScan | Pfam, TIGRFAM, Gene3D, PANTHER, CDD, SUPERFAMILY, PRINTS | BLAST returns LOW or hypothetical |
| 3 | RCSB PDB search | PDB experimental + AlphaFold CSMs | BLAST and InterProScan both fail |
| 4 | Composition analysis | Deterministic, no API | All three tiers return LOW or UNKNOWN, sequence ≤400 aa |

### Confidence tiers

- `HIGH` — strong homology to a named, characterized protein
- `MODERATE` — meaningful signal but incomplete (hypothetical hit, partial coverage)
- `LOW` — weak evidence or composition-only hypothesis
- `UNKNOWN` — no hits above threshold from this tool

### Dark matter scenarios

When all three database tiers fail, ORACLE distinguishes two biologically different outcomes:

**Scenario A** — BLAST found a homolog but it is also uncharacterized:
> `Uncharacterized protein (homolog detected in Vicingus sp.)`

This means the protein is conserved across at least two genomes, implying selective pressure and likely function. It is dark matter with relatives.

**Scenario B** — Nothing in any database:
> `Uncharacterized protein (no database representatives found)`

This is a genuine orphan with no sequence, domain, or structural relatives. It may represent a novel fold or lineage-specific innovation.

## Example results

Two sequences from the BATS (Bermuda Atlantic Time-series Study) metagenomic dataset demonstrate a piece of the uninformative and uncharacterized annotation spectrum. 

### Weak hit — BATS hypothetical protein

```
Sequence:   Ga0531459_000429_2805_3143_5  (112 aa)
Result:     Uncharacterized protein (homolog detected in Vicingus sp.)
Confidence: LOW
Tools used: BLAST → InterProScan → RCSB PDB → Composition
```

BLAST finds a 90% coverage hit to a hypothetical protein in *Vicingus* phage (e-value 2.1e-5). InterProScan detects only intrinsic disorder (47% MobiDB-lite). RCSB finds no structural relatives. Composition analysis identifies a basic N-terminus (30% K/R/H, consistent with DNA binding) and an acidic C-terminus (25% D/E, consistent with protein interaction).

### Dark matter — BATS orphan protein

```
Sequence:   Ga0531459_000428_12414_11821_6  (197 aa)
Result:     Uncharacterized protein (no database representatives found)
Confidence: LOW
Tools used: BLAST → InterProScan → RCSB PDB → Composition
```

No BLAST hits. InterProScan detects only a weak coiled-coil signal (14% coverage, e-value 1.0). RCSB finds nothing in experimental PDB or AlphaFold CSMs. Composition analysis finds high disorder propensity (disorder index 0.41) with enrichment in disorder-promoting residues, consistent with a regulatory IDP or hub protein.

### Animated protein cards

For LOW or UNKNOWN confidence sequences, ORACLE generates an animated HTML card 
with a spinning ESMFold predicted structure colored by pLDDT confidence.

| Sequence | Result | Card |
|----------|--------|------|
| Ga0531459_000428 (197 aa) | Uncharacterized — no database representatives | [View card ↗](https://zschreib.github.io/oracle/example_output/protein_cards/Ga0531459_000428_12414_11821_6_20260525_184209.card.html) |
| Ga0531459_000429 (112 aa) | Uncharacterized — homolog in *Vicingus sp.* | [View card ↗](https://zschreib.github.io/oracle/example_output/protein_cards/Ga0531459_000429_2805_3143_5_20260525_182650.card.html) |

---

## Installation

**Conda (recommended for bioinformatics environments):**

```bash
git clone https://github.com/zschreib/oracle.git
cd oracle
conda env create -f environment.yml
conda activate oracle
pip install -e .
```

**pip:**

```bash
git clone https://github.com/zschreib/oracle.git
cd oracle
pip install -e .
```

Copy the environment template and fill in your keys:

```bash
cp .env.example .env
```

```bash
# .env
# Anthropic API key — create at: https://console.anthropic.com/settings/keys
ANTHROPIC_API_KEY=

# NCBI API key — raises BLAST rate limit from 3 to 10 requests/second
# Register free at: https://www.ncbi.nlm.nih.gov/account/
NCBI_API_KEY=
```

---
## Usage

```bash
# Annotate a single sequence
oracle -i sequence.fasta

# Verbose output showing agent reasoning at each step
oracle -i sequence.fasta --verbose

# Dry run (shows tool schemas, no API calls)
oracle -i sequence.fasta --dry-run
```

Output is written to `reports/` as both `.report.txt` and `.report.json`. 

### Generate an animated protein card

For LOW or UNKNOWN confidence sequences under 400 aa, generate an animated HTML card with a spinning ESMFold predicted structure colored by pLDDT:

```bash
python generate_card.py reports/output.report.json --fasta example_input/input_sequence.fasta
```

The card is self-contained HTML with no runtime dependencies. Open in any browser.

> The card is intended for dark matter sequences where ESMFold structure prediction is the only structural hypothesis available. For well-characterized sequences, the annotation report is the primary output.

## Project structure

```
oracle/
  oracle/
    agent/
      oracle.py       ← agent loop, tool dispatch, Anthropic tool use API
      prompts.py      ← system prompt and escalation rules
    tools/
      blast.py        ← BLASTp via NCBI REST API (tier 1)
      hmmer.py        ← InterProScan via EBI Job Dispatcher (tier 2)
      structure.py    ← RCSB PDB search + AlphaFold CSM fallback (tier 3)
      composition.py  ← deterministic sequence composition analysis (tier 4)
      registry.py     ← tool registration and dispatch
      base.py         ← AnnotationTool abstract base class
    utils/
      report.py       ← render OracleReport to text and JSON
    main.py           ← CLI entry point
    models.py         ← ConfidenceTier, AnnotationEvidence, AgentStep, OracleReport
  tests/
    mock.py           ← MockTool and test scenarios
    test_rcsb_api.py  ← RCSB PDB API tests including CSM fallback
    test_hmmer_api.py ← EBI InterProScan API tests
    test_esmfold_api.py ← ESMFold structure prediction tests
  example_input/      
  generate_card.py    ← animated HTML card generator for dark matter sequences

```

---

## Composition analysis rules

The tier 4 composition tool fires only when all three database tiers return LOW or UNKNOWN confidence and the sequence is under 400 aa. It applies deterministic biophysical rules with three key edge case guards:

**Signal peptide guard:** A basic N-terminus (≥25% K/R/H in first 30 aa) followed by a hydrophobic h-region (≥8 aa run in positions 6-25) is flagged as a signal peptide rather than a DNA-binding domain. The downstream hydrophobic stretch is not double-counted as a transmembrane helix.

**Transmembrane threshold:** A hydrophobic run of 15+ aa is required for a transmembrane helix call. Runs of 10-14 aa are noted as potential hydrophobic core of a soluble protein without making a membrane protein call.

**Cysteine context:** Two cysteines within 30 aa are flagged as `disulfide_or_metal` with explicit text noting that interpretation depends on environment: oxidizing (phage virion, periplasm) favors disulfide; reducing (cytoplasm, strict anaerobes) favors metal coordination. Three or more cysteines are flagged as metal-binding since that is the more likely interpretation regardless of context.

---

## Requirements

- Python 3.11+
- Anthropic API key
- NCBI API key (free, increases BLAST rate limit)
- `anthropic>=0.25.0`, `requests>=2.31.0`, `biopython>=1.83`

BLAST and InterProScan run against remote APIs. No local database installation required.

---

## Summary

ORACLE was built to address a specific gap in metagenomic annotation: existing pipelines assign "hypothetical protein" to dark matter sequences and stop there. The composition analysis tier in particular was motivated by phage dark matter work from the BATS dataset, where short proteins with bipartite charge distributions are systematically missed by domain-based annotation approaches.

The two example sequences are real proteins from the BATS dataset used in metagenomic studies of oceanic viral communities.

---

## Disclaimer

ORACLE is a research tool intended to assist in hypothesis generation for uncharacterized protein sequences. All annotations, especially those at LOW or UNKNOWN confidence, are computational predictions and should not be treated as experimentally validated functional assignments.

Composition-based annotations (tier 4) are speculative by design. They are derived from physicochemical properties of the amino acid sequence and have no direct experimental support. They are provided to narrow the hypothesis space for follow-up experiments, not to replace them.

ORACLE makes live calls to NCBI BLAST, EBI InterProScan, and RCSB PDB. Results may vary over time as these databases are updated. Annotations are not guaranteed to be reproducible across database versions.

This agent is not intended for clinical, diagnostic, or regulatory use.


## Author

zschreib.dev@gmail.com