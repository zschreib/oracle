"""
Sequence composition analysis tool for ORACLE.

A deterministic, rule-based analysis of protein sequence properties
that fires only when BLAST, InterProScan, and RCSB PDB all return
LOW or UNKNOWN confidence. Produces structured functional hypotheses
from sequence composition rather than database comparison.

This tool never claims HIGH confidence because composition alone
cannot confirm function. Its value is in narrowing the hypothesis
space for dark matter sequences and providing reproducible,
interpretable evidence for experimental follow-up.

Rules implemented here are based on well-established biophysical
principles and are directly relevant to phage and metagenomic
dark matter proteins:

- Basic N-terminus: K/R enrichment in first 30aa suggests DNA/RNA binding
- Acidic C-terminus: D/E enrichment in last 20aa suggests protein interaction
- Hydrophobic core: extended hydrophobic stretches suggest transmembrane or
  signal peptides
- Cysteine clustering: multiple cysteines within 30aa suggests disulfide
  bonds or metal coordination
- Repeat structure: tandem repeats suggest structural or scaffolding roles
- Charge distribution: overall pI predicts cellular localization tendencies
- Small protein bias: <100aa proteins in phage contexts are often
  regulatory, anti-restriction, or virion structural accessories

No external API calls.
"""

from .base import AnnotationTool
from ..models import AnnotationEvidence, ConfidenceTier


# Amino acid property sets
_BASIC     = set("KRH")
_ACIDIC    = set("DE")
_HYDROPHOB = set("VILMFYW")
_POLAR     = set("STCNQ")
_TINY      = set("AGSP")


class CompositionTool(AnnotationTool):
    """
    Rule-based sequence composition analysis for dark matter proteins.

    Only called when all upstream tools return LOW or UNKNOWN.
    Returns LOW confidence functional hypotheses based on sequence
    properties. Never returns HIGH or MODERATE confidence since
    composition alone cannot confirm function.
    """

    @property
    def name(self) -> str:
        return "composition"

    @property
    def tool_schema(self) -> dict:
        return {
            "name": self.name,
            "description": (
                "Run deterministic sequence composition analysis to generate "
                "functional hypotheses for uncharacterized proteins. "
                "ONLY use this when blast_nr, hmmer_pfam, AND rcsb_pdb have all "
                "returned LOW or UNKNOWN confidence. "
                "This tool never returns HIGH or MODERATE confidence. "
                "It provides structured interpretable evidence for experimental "
                "follow-up on genuine dark matter sequences. "
                "Runs instantly, no API calls."
            ),
            "input_schema": {
                "type": "object",
                "properties": {
                    "sequence": {
                        "type": "string",
                        "description": "Amino acid sequence string, no header, no whitespace."
                    }
                },
                "required": ["sequence"]
            }
        }

    def run(self, sequence: str, **kwargs) -> AnnotationEvidence:
        """
        Analyze sequence composition and return structured hypotheses.

        Runs all composition rules and produces a ranked list of
        functional hypotheses with supporting evidence.
        """
        sequence = sequence.strip().replace("\n", "").replace(" ", "").upper()

        if len(sequence) < 10:
            return self._failure("Sequence too short for composition analysis")

        findings = []
        hypotheses = []

        # Run all composition rules
        findings += self._analyze_termini(sequence)
        findings += self._analyze_hydrophobicity(sequence, prior_findings=findings)
        findings += self._analyze_cysteines(sequence)
        findings += self._analyze_repeats(sequence)
        findings += self._analyze_size_context(sequence)
        findings += self._analyze_charge(sequence)
        findings += self._analyze_disorder_proxies(sequence)

        # Build hypothesis list from findings
        for f in findings:
            if f.get("hypothesis"):
                hypotheses.append(f["hypothesis"])

        # Build annotation string
        if hypotheses:
            annotation = "Uncharacterized protein"
            if len(hypotheses) == 1:
                annotation = f"Putative {hypotheses[0]}"
            else:
                annotation = f"Putative {hypotheses[0]} (alternative: {hypotheses[1]})"
        else:
            annotation = "Uncharacterized protein, no compositional signal"

        # Build reasoning string from all findings
        reasoning_parts = []
        for f in findings:
            if f.get("detail"):
                reasoning_parts.append(f["detail"])

        reasoning = (
            f"Sequence composition analysis ({len(sequence)} aa). "
            + " ".join(reasoning_parts)
            + " All annotations are speculative and require experimental validation."
        )

        return AnnotationEvidence(
            tool_name=self.name,
            annotation=annotation,
            confidence=ConfidenceTier.LOW,
            score=None,
            coverage=None,
            hit_id=None,
            hit_description="; ".join(hypotheses) if hypotheses else "no signal",
            reasoning=reasoning,
            raw_output=str(findings)
        )

    def _analyze_termini(self, seq: str) -> list:
        """
        Check N and C terminal charge enrichment.

        Basic N-terminus: DNA/RNA binding, nuclear localization, membrane
        translocation. Acidic C-terminus: protein-protein interaction hub,
        regulatory tail, intrinsic disorder anchor.

        Signal peptide guard: a basic N-terminus (n-region) immediately
        followed by a hydrophobic stretch (h-region) is the classic
        Sec-dependent signal peptide pattern. In that case we suppress
        the DNA-binding hypothesis and flag secretion instead.
        """
        findings = []
        n_region = seq[:30]
        c_region = seq[-20:]

        # N-terminal basic residues
        basic_n = sum(1 for aa in n_region if aa in _BASIC)
        basic_density = basic_n / len(n_region)

        if basic_density >= 0.25:
            # Signal peptide guard: check for hydrophobic h-region in
            # positions 6-25 (after the charged n-region)
            h_region = seq[5:25]
            h_run = 0
            max_h_run = 0
            for aa in h_region:
                if aa in _HYDROPHOB:
                    h_run += 1
                    max_h_run = max(max_h_run, h_run)
                else:
                    h_run = 0

            if max_h_run >= 8:
                # Looks like a signal peptide, not a DNA-binding domain
                findings.append({
                    "rule": "signal_peptide",
                    "value": max_h_run,
                    "hypothesis": "secreted or periplasmic protein",
                    "detail": (
                        f"Signal peptide candidate: basic n-region "
                        f"({basic_density:.0%} K/R/H in first 30aa) followed by "
                        f"hydrophobic h-region (run of {max_h_run} in positions 6-25). "
                        f"Classic Sec-dependent signal peptide pattern. "
                        f"Protein may be secreted or periplasmic rather than "
                        f"cytoplasmic DNA-binding."
                    )
                })
            else:
                # No h-region, genuine basic N-terminus
                findings.append({
                    "rule": "basic_n_terminus",
                    "value": round(basic_density, 2),
                    "hypothesis": "DNA/RNA-binding protein",
                    "detail": (
                        f"Basic N-terminus: {basic_density:.0%} K/R/H density "
                        f"in first 30aa (>{0.25:.0%} threshold), no downstream "
                        f"hydrophobic h-region detected. "
                        f"Consistent with nucleic acid binding, membrane translocation, "
                        f"or phage DNA-injection machinery."
                    )
                })

        # C-terminal acidic residues
        acidic_c = sum(1 for aa in c_region if aa in _ACIDIC)
        acidic_density = acidic_c / len(c_region)

        if acidic_density >= 0.25:
            findings.append({
                "rule": "acidic_c_terminus",
                "value": round(acidic_density, 2),
                "hypothesis": "protein-protein interaction mediator",
                "detail": (
                    f"Acidic C-terminus: {acidic_density:.0%} D/E density "
                    f"in last 20aa. "
                    f"Consistent with intrinsically disordered interaction tail, "
                    f"acidic activation domain, or hub protein interface."
                )
            })

        return findings

    def _analyze_hydrophobicity(self, seq: str, prior_findings: list = None) -> list:
        """
        Detect extended hydrophobic stretches suggesting membrane association.

        Transmembrane helices require 15-25 hydrophobic residues to span
        a lipid bilayer (~30 Angstroms at 1.5 Å per residue in helix).
        Shorter hydrophobic stretches can be interior of soluble globular
        proteins and should not be flagged as transmembrane.

        Signal peptide guard: if a signal peptide was already detected
        in the N-terminal region, a hydrophobic stretch there is the
        expected h-region and should not be double-counted as a
        transmembrane helix.
        """
        findings = []
        prior_findings = prior_findings or []
        signal_peptide_flagged = any(
            f.get("rule") == "signal_peptide" for f in prior_findings
        )

        # Find all hydrophobic runs and their positions
        runs = []
        current_run = 0
        current_start = 0
        for i, aa in enumerate(seq):
            if aa in _HYDROPHOB:
                if current_run == 0:
                    current_start = i
                current_run += 1
            else:
                if current_run > 0:
                    runs.append((current_start, current_start + current_run, current_run))
                current_run = 0
        if current_run > 0:
            runs.append((current_start, current_start + current_run, current_run))

        if not runs:
            return findings

        max_run_start, max_run_end, max_run = max(runs, key=lambda r: r[2])

        # Suppress N-terminal hydrophobic runs if signal peptide already flagged
        if signal_peptide_flagged and max_run_start < 30:
            # Check if any non-N-terminal run is long enough
            other_runs = [(s, e, l) for s, e, l in runs if s >= 30 and l >= 15]
            if not other_runs:
                return findings
            max_run_start, max_run_end, max_run = max(other_runs, key=lambda r: r[2])

        # 15aa minimum for transmembrane helix (bilayer crossing requirement)
        if max_run >= 15:
            findings.append({
                "rule": "transmembrane",
                "value": max_run,
                "hypothesis": "membrane protein or phage spanin",
                "detail": (
                    f"Transmembrane candidate: hydrophobic run of {max_run} residues "
                    f"at position {max_run_start+1}-{max_run_end} "
                    f"(≥15aa required to span lipid bilayer). "
                    f"Consistent with integral membrane protein, phage spanin, "
                    f"or holin family protein."
                )
            })
        elif max_run >= 10:
            # Flag but note it could be soluble protein interior
            findings.append({
                "rule": "hydrophobic_core",
                "value": max_run,
                "hypothesis": None,
                "detail": (
                    f"Moderate hydrophobic run of {max_run} residues at "
                    f"position {max_run_start+1}-{max_run_end}. "
                    f"Below the 15aa transmembrane threshold. May be hydrophobic "
                    f"core of a soluble globular protein, amphipathic helix, "
                    f"or short signal anchor."
                )
            })

        return findings

    def _analyze_cysteines(self, seq: str) -> list:
        """
        Detect cysteine clustering suggesting disulfide bonds or metal binding.

        Interpretation depends on cellular context which is unknown for dark
        matter sequences, so both hypotheses are reported:
        - Extracellular/phage virion context: disulfide bonds (oxidizing)
        - Intracellular/anaerobic context: metal coordination (reducing)

        Three or more cysteines with clustering strongly suggests metal
        coordination (zinc finger, iron-sulfur cluster) over simple disulfide,
        since disulfide bonds rarely involve more than two cysteines per domain.
        """
        findings = []
        cys_positions = [i for i, aa in enumerate(seq) if aa == "C"]

        if len(cys_positions) < 2:
            return findings

        clusters = []
        for i in range(len(cys_positions) - 1):
            if cys_positions[i+1] - cys_positions[i] <= 30:
                clusters.append((cys_positions[i], cys_positions[i+1]))

        total_cys = len(cys_positions)

        if total_cys >= 3 and clusters:
            findings.append({
                "rule": "metal_binding",
                "value": total_cys,
                "hypothesis": "metal-binding or redox-active protein",
                "detail": (
                    f"Cysteine cluster: {total_cys} cysteines total, "
                    f"{len(clusters)} pair(s) within 30aa. "
                    f"Three or more clustered cysteines strongly suggests metal "
                    f"coordination (zinc finger, iron-sulfur cluster) in aerobic "
                    f"contexts, or active-site cysteine in redox enzymes in anaerobic "
                    f"contexts. Disulfide bonds less likely with this many cysteines."
                )
            })
        elif total_cys >= 2 and clusters:
            findings.append({
                "rule": "disulfide_or_metal",
                "value": total_cys,
                "hypothesis": "disulfide-stabilized or metal-coordinating protein",
                "detail": (
                    f"Cysteine pair within 30aa (positions "
                    f"{clusters[0][0]+1} and {clusters[0][1]+1}). "
                    f"In oxidizing environments (phage virion exterior, periplasm): "
                    f"likely disulfide bond for structural stability. "
                    f"In reducing environments (cytoplasm, strict anaerobes): "
                    f"likely metal coordination or catalytic cysteine. "
                    f"Context unknown for this dark matter sequence."
                )
            })

        return findings

    def _analyze_repeats(self, seq: str) -> list:
        """
        Detect simple tandem repeats suggesting structural or scaffolding roles.

        Looks for repeated subsequences of length 3-8 appearing 3+ times.
        Common in structural phage proteins (tail fibers, baseplate components).
        """
        findings = []

        for repeat_len in range(4, 9):
            counts = {}
            for i in range(len(seq) - repeat_len):
                unit = seq[i:i+repeat_len]
                # Skip low-complexity units (all same AA)
                if len(set(unit)) < 2:
                    continue
                counts[unit] = counts.get(unit, 0) + 1

            top_repeats = [(u, c) for u, c in counts.items() if c >= 3]
            if top_repeats:
                best = max(top_repeats, key=lambda x: x[1])
                findings.append({
                    "rule": "tandem_repeat",
                    "value": best[1],
                    "hypothesis": "repeat-containing structural protein",
                    "detail": (
                        f"Tandem repeat detected: '{best[0]}' appears {best[1]} times. "
                        f"Consistent with structural scaffold, tail fiber protein, "
                        f"or repeat-containing binding protein."
                    )
                })
                break  # One repeat finding is enough

        return findings

    def _analyze_size_context(self, seq: str) -> list:
        """
        Use protein length to constrain functional hypotheses.

        Very small proteins (<80aa) in phage contexts have characteristic
        functional categories. Size combined with other signals narrows
        the hypothesis space significantly.
        """
        findings = []
        n = len(seq)

        if n < 80:
            findings.append({
                "rule": "very_small",
                "value": n,
                "hypothesis": None,
                "detail": (
                    f"Very small protein ({n} aa). "
                    f"Phage proteins of this size are often anti-restriction "
                    f"factors, superinfection exclusion proteins, DNA injection "
                    f"pilots, or virion structural accessories."
                )
            })
        elif n < 150:
            findings.append({
                "rule": "small",
                "value": n,
                "hypothesis": None,
                "detail": (
                    f"Small protein ({n} aa). "
                    f"Consistent with regulatory protein, single-domain enzyme, "
                    f"or structural accessory."
                )
            })

        return findings

    def _analyze_charge(self, seq: str) -> list:
        """
        Estimate isoelectric point and overall charge character.

        Highly basic proteins (pI > 9) often interact with nucleic acids.
        Highly acidic proteins (pI < 5) are often intrinsically disordered
        or involved in signaling/regulatory roles.
        """
        findings = []
        pos = sum(1 for aa in seq if aa in "KRH")
        neg = sum(1 for aa in seq if aa in "DE")
        net_charge = pos - neg
        charge_ratio = net_charge / len(seq)

        if charge_ratio > 0.15:
            findings.append({
                "rule": "highly_basic",
                "value": round(charge_ratio, 2),
                "hypothesis": None,
                "detail": (
                    f"Highly basic overall charge (net charge ratio {charge_ratio:+.2f}). "
                    f"Predicted pI > 9. Consistent with nucleic acid-binding, "
                    f"ribosome-associated, or membrane-active protein."
                )
            })
        elif charge_ratio < -0.15:
            findings.append({
                "rule": "highly_acidic",
                "value": round(charge_ratio, 2),
                "hypothesis": None,
                "detail": (
                    f"Highly acidic overall charge (net charge ratio {charge_ratio:+.2f}). "
                    f"Predicted pI < 5. Consistent with intrinsically disordered "
                    f"regulatory protein or acidic activation domain."
                )
            })

        return findings

    def _analyze_disorder_proxies(self, seq: str) -> list:
        """
        Use compositional proxies for intrinsic disorder.

        Disorder-promoting residues: E, K, R, S, Q, P, A, G (PEST regions).
        Order-promoting residues: W, Y, F, I, L, V, C, M, H, N.
        High disorder-promoting fraction suggests IDP characteristics.
        """
        findings = []
        disorder_promoting = set("EKRSQPAG")
        order_promoting = set("WYFILVCMHN")

        dp = sum(1 for aa in seq if aa in disorder_promoting) / len(seq)
        op = sum(1 for aa in seq if aa in order_promoting) / len(seq)
        disorder_index = dp - op

        if disorder_index > 0.25:
            findings.append({
                "rule": "high_disorder",
                "value": round(disorder_index, 2),
                "hypothesis": "intrinsically disordered regulatory protein",
                "detail": (
                    f"High disorder propensity (disorder index {disorder_index:.2f}). "
                    f"Enriched in disorder-promoting residues (E/K/R/S/Q/P/A/G). "
                    f"Consistent with hub protein, regulatory IDP, or flexible linker."
                )
            })

        return findings

    def _failure(self, message: str) -> AnnotationEvidence:
        return AnnotationEvidence(
            tool_name=self.name,
            annotation="Tool failure",
            confidence=ConfidenceTier.UNKNOWN,
            reasoning=message
        )