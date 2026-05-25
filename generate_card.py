"""
ORACLE protein card generator.

Generates an animated HTML card for uncharacterized or low-confidence
sequences where visual structure prediction adds scientific value.

INTENDED USE CASE
-----------------
This card is designed for sequences that ORACLE annotates with LOW or
UNKNOWN confidence, typically dark matter proteins with no database
representatives. For these sequences, an ESMFold predicted structure
provides the only structural hypothesis available and is scientifically
meaningful to show.

For HIGH confidence sequences (strong BLAST hits, known proteins), the
annotation report itself is the primary output. Those proteins already
have experimental structures in PDB and a predicted structure adds little.

SEQUENCE LENGTH LIMIT
---------------------
ESMFold API accepts sequences up to ~400 amino acids. The card generator
will exit with an informative message if the sequence exceeds this limit.
Use the annotation report for longer well-characterized sequences.

Usage:
    python generate_card.py reports/seq_001.report.json --fasta example_input/dark_matter.fasta
    python generate_card.py reports/seq_001.report.json --sequence MKTII...

Output:
    reports/seq_001.card.html  — self-contained animated HTML, no runtime dependencies
"""

import argparse
import json
import os
import re
import sys
import requests


ESMFOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"

# pLDDT color scheme matching AlphaFold's official visualization
# Very high (90+): blue, High (70-90): cyan, Medium (50-70): yellow, Low (<50): orange
PLDDT_COLORS = [
    (0,   "#ff7d45"),   # <50  orange
    (50,  "#ffdb13"),   # 50-70 yellow
    (70,  "#65cbf3"),   # 70-90 cyan
    (90,  "#0053d6"),   # 90+  blue
]


def load_report(report_path: str) -> dict:
    with open(report_path) as f:
        return json.load(f)


def load_sequence(fasta_path: str) -> tuple[str, str]:
    """Return (sequence_id, sequence) from FASTA file."""
    seq_id = ""
    seq = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                seq_id = line[1:].split()[0]
            else:
                seq.append(line.upper())
    return seq_id, "".join(seq)


def fold_sequence(sequence: str) -> str:
    """Call ESMFold API, return PDB string."""
    print(f"  Predicting structure ({len(sequence)} aa)...")
    response = requests.post(
        ESMFOLD_URL,
        data=sequence,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        timeout=60
    )
    response.raise_for_status()
    return response.text


def parse_pdb_for_js(pdb_str: str) -> str:
    """
    Extract Cα atoms from PDB and format as a JS array.

    Each entry: [x, y, z, plddt] where plddt is 0-100 from B-factor column.
    Returns a JSON-serializable string for embedding in HTML.
    """
    atoms = []
    for line in pdb_str.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                plddt = float(line[60:66])  # B-factor = pLDDT in ESMFold
                atoms.append([round(x, 2), round(y, 2), round(z, 2), round(plddt * 100 if plddt <= 1.0 else plddt, 1)])
            except (ValueError, IndexError):
                continue
    return json.dumps(atoms)


def compute_plddt_stats(pdb_str: str) -> dict:
    """Compute mean pLDDT and per-residue counts by confidence tier."""
    plddts = []
    for line in pdb_str.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                plddts.append(float(line[60:66]))
            except (ValueError, IndexError):
                continue

    if not plddts:
        return {"mean": 0, "very_high": 0, "high": 0, "medium": 0, "low": 0, "total": 0}

    # ESMFold reports pLDDT in 0-100 range in the B-factor column.
    # Normalize to 0-100 if values appear to be in 0-1 range.
    if max(plddts) <= 1.0:
        plddts = [p * 100 for p in plddts]

    mean = sum(plddts) / len(plddts)
    return {
        "mean": round(mean, 1),
        "very_high": sum(1 for p in plddts if p >= 90),
        "high":      sum(1 for p in plddts if 70 <= p < 90),
        "medium":    sum(1 for p in plddts if 50 <= p < 70),
        "low":       sum(1 for p in plddts if p < 50),
        "total":     len(plddts)
    }


def confidence_color(tier: str) -> str:
    return {
        "HIGH":     "#4ade80",
        "MODERATE": "#facc15",
        "LOW":      "#fb923c",
        "UNKNOWN":  "#94a3b8",
    }.get(tier, "#94a3b8")


def confidence_bar(tier: str) -> str:
    fills = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "UNKNOWN": 1}
    n = fills.get(tier, 1)
    color = confidence_color(tier)
    bars = "".join(
        f'<span style="background:{color if i < n else "rgba(255,255,255,0.1)"}">'
        f'</span>'
        for i in range(4)
    )
    return f'<div class="conf-bars">{bars}</div>'


def generate_html(report: dict, pdb_str: str, plddt: dict) -> str:
    """Generate the full self-contained animated HTML card."""

    atoms_json = parse_pdb_for_js(pdb_str)
    tier = report.get("final_confidence", "UNKNOWN")
    annotation = report.get("final_annotation", "Uncharacterized protein")
    seq_id = report.get("sequence_id", "")
    seq_len = report.get("sequence_length", 0)
    warnings = report.get("warnings", [])

    # Build evidence chain summary rows
    evidence_rows = ""
    for step in report.get("evidence_chain", []):
        tool = step.get("tool", "").upper()
        conf = step.get("confidence", "UNKNOWN")
        ann = step.get("annotation", "")
        score = step.get("score")
        if score is None:
            score_str = "—"
        elif score == 0.0:
            score_str = "0.0"
        elif score >= 1.0:
            score_str = "n/s"   # not significant
        elif score < 0.001:
            score_str = f"{score:.1e}"
        else:
            score_str = f"{score:.3f}"
        col = confidence_color(conf)
        evidence_rows += f"""
        <tr>
          <td class="tool-name">{tool}</td>
          <td class="tool-ann">{ann}</td>
          <td class="tool-tier"><span class="badge" style="background:{col}20;color:{col};border:1px solid {col}40">{conf}</span></td>
          <td class="score-val">{score_str}</td>
        </tr>"""

    skipped_rows = ""
    for skip in report.get("skipped_tools", []):
        tool = skip.get("tool", "").upper()
        skipped_rows += f"""
        <tr class="skipped">
          <td class="tool-name">{tool}</td>
          <td class="tool-ann" colspan="3">skipped — evidence sufficient</td>
        </tr>"""

    plddt_bar_html = ""
    if plddt["total"] > 0:
        for label, key, color in [
            ("≥90", "very_high", "#0053d6"),
            ("70–90", "high", "#65cbf3"),
            ("50–70", "medium", "#ffdb13"),
            ("<50", "low", "#ff7d45"),
        ]:
            pct = plddt[key] / plddt["total"] * 100
            if pct > 0:
                plddt_bar_html += (
                    f'<div style="width:{pct:.1f}%;background:{color};'
                    f'height:100%;display:inline-block;"></div>'
                )

    warning_html = ""
    if warnings:
        warning_html = '<div class="warnings">' + "".join(
            f'<div class="warn-item">⚠ {w}</div>' for w in warnings
        ) + "</div>"

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>ORACLE — {seq_id}</title>
<style>
  @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@300;400;500;600&display=swap');

  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}

  body {{
    background: #080c14;
    min-height: 100vh;
    display: flex;
    align-items: center;
    justify-content: center;
    font-family: 'IBM Plex Sans', sans-serif;
    padding: 2rem;
  }}

  .card {{
    width: 1100px;
    max-width: 100%;
    background: #0d1520;
    border: 1px solid #1e2d42;
    border-radius: 16px;
    display: grid;
    grid-template-columns: 1fr 460px;
    overflow: hidden;
    box-shadow: 0 0 80px rgba(0,120,255,0.08), 0 32px 64px rgba(0,0,0,0.6);
  }}

  /* LEFT PANEL */
  .left {{
    padding: 1.75rem 2rem;
    border-right: 1px solid #1e2d42;
    display: flex;
    flex-direction: column;
    gap: 1rem;
  }}

  .oracle-badge {{
    display: inline-flex;
    align-items: center;
    gap: 0.4rem;
    font-family: 'IBM Plex Mono', monospace;
    font-size: 10px;
    font-weight: 500;
    letter-spacing: 0.12em;
    color: #4a7fa5;
    text-transform: uppercase;
  }}

  .oracle-badge::before {{
    content: '';
    width: 6px; height: 6px;
    border-radius: 50%;
    background: #4a7fa5;
    animation: pulse 2s ease-in-out infinite;
  }}

  @keyframes pulse {{
    0%, 100% {{ opacity: 1; transform: scale(1); }}
    50% {{ opacity: 0.4; transform: scale(0.8); }}
  }}

  .seq-id {{
    font-family: 'IBM Plex Mono', monospace;
    font-size: 11px;
    color: #3d5a73;
    word-break: break-all;
  }}

  .annotation {{
    font-size: 1.1rem;
    font-weight: 600;
    color: #e2eaf4;
    line-height: 1.35;
    letter-spacing: -0.01em;
  }}

  .conf-row {{
    display: flex;
    align-items: center;
    gap: 0.75rem;
  }}

  .conf-label {{
    font-family: 'IBM Plex Mono', monospace;
    font-size: 11px;
    color: #4a7fa5;
    letter-spacing: 0.08em;
    text-transform: uppercase;
  }}

  .conf-bars {{
    display: flex;
    gap: 3px;
  }}

  .conf-bars span {{
    width: 18px;
    height: 6px;
    border-radius: 2px;
  }}

  .conf-tier {{
    font-family: 'IBM Plex Mono', monospace;
    font-size: 11px;
    font-weight: 500;
    letter-spacing: 0.06em;
  }}

  .divider {{
    height: 1px;
    background: #1e2d42;
  }}

  .section-label {{
    font-family: 'IBM Plex Mono', monospace;
    font-size: 10px;
    color: #3d5a73;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    margin-bottom: 0.5rem;
  }}

  table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 12px;
    table-layout: fixed;
  }}

  td {{
    padding: 5px 4px;
    color: #8aa4bc;
    vertical-align: top;
  }}

  .tool-name {{
    font-family: 'IBM Plex Mono', monospace;
    font-size: 9px;
    color: #4a7fa5;
    width: 80px;
    letter-spacing: 0.03em;
    white-space: nowrap;
  }}

  .tool-ann {{
    color: #c8d8e8;
    font-size: 11px;
    word-break: break-word;
    line-height: 1.4;
  }}

  .tool-tier {{
    width: 72px;
    white-space: nowrap;
  }}

  .score-val {{
    font-family: 'IBM Plex Mono', monospace;
    font-size: 10px;
    color: #3d5a73;
    text-align: right;
    width: 70px;
    white-space: nowrap;
  }}

  .badge {{
    font-family: 'IBM Plex Mono', monospace;
    font-size: 9px;
    padding: 2px 6px;
    border-radius: 4px;
    white-space: nowrap;
    letter-spacing: 0.05em;
  }}

  .conf-legend {{
    display: flex;
    flex-wrap: wrap;
    gap: 0.3rem 0.6rem;
    margin-top: 0.5rem;
  }}

  .conf-legend > span {{
    font-size: 9px;
    color: #3d5a73;
    display: flex;
    align-items: center;
    gap: 4px;
    font-family: 'IBM Plex Mono', monospace;
  }}

  .conf-legend .badge {{
    font-size: 8px;
    padding: 1px 4px;
    border-radius: 3px;
    font-weight: 500;
    letter-spacing: 0.04em;
    flex-shrink: 0;
  }}

  tr.skipped td {{
    opacity: 0.35;
    font-style: italic;
  }}

  .plddt-section {{ display: flex; flex-direction: column; gap: 0.4rem; }}

  .plddt-bar-wrap {{
    height: 6px;
    border-radius: 3px;
    overflow: hidden;
    background: rgba(255,255,255,0.05);
  }}

  .plddt-mean {{
    font-family: 'IBM Plex Mono', monospace;
    font-size: 11px;
    color: #65cbf3;
  }}

  .plddt-legend {{
    display: flex;
    gap: 0.75rem;
    flex-wrap: wrap;
  }}

  .plddt-legend span {{
    font-size: 10px;
    color: #4a7fa5;
    display: flex;
    align-items: center;
    gap: 3px;
  }}

  .plddt-legend span::before {{
    content: '';
    width: 8px; height: 8px;
    border-radius: 2px;
    display: inline-block;
  }}

  .warnings {{
    display: flex;
    flex-direction: column;
    gap: 0.3rem;
  }}

  .warn-item {{
    font-size: 10px;
    color: #b45309;
    background: rgba(180, 83, 9, 0.08);
    border: 1px solid rgba(180, 83, 9, 0.2);
    border-radius: 4px;
    padding: 4px 8px;
    line-height: 1.4;
  }}

  .meta-row {{
    display: flex;
    gap: 1rem;
    flex-wrap: wrap;
  }}

  .meta-item {{
    font-family: 'IBM Plex Mono', monospace;
    font-size: 10px;
    color: #3d5a73;
  }}

  .meta-item strong {{
    color: #4a7fa5;
    font-weight: 500;
  }}

  /* RIGHT PANEL */
  .right {{
    position: relative;
    background: #080c14;
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    padding: 1rem;
    min-height: 520px;
  }}

  canvas {{
    width: 100% !important;
    height: 100% !important;
    border-radius: 8px;
  }}

  .structure-label {{
    position: absolute;
    bottom: 1rem;
    left: 50%;
    transform: translateX(-50%);
    font-family: 'IBM Plex Mono', monospace;
    font-size: 9px;
    color: #1e2d42;
    letter-spacing: 0.1em;
    text-transform: uppercase;
    white-space: nowrap;
  }}

  .plddt-legend-3d {{
    position: absolute;
    top: 1rem;
    right: 1rem;
    display: flex;
    flex-direction: column;
    gap: 3px;
  }}

  .plddt-legend-3d span {{
    font-size: 9px;
    font-family: 'IBM Plex Mono', monospace;
    color: #3d5a73;
    display: flex;
    align-items: center;
    gap: 4px;
  }}

  .plddt-legend-3d span::before {{
    content: '';
    width: 8px; height: 3px;
    border-radius: 2px;
    display: inline-block;
  }}
</style>
</head>
<body>
<div class="card">

  <!-- LEFT -->
  <div class="left">
    <div>
      <div class="oracle-badge">ORACLE Annotation Report</div>
      <div class="seq-id" style="margin-top:0.4rem">{seq_id}</div>
    </div>

    <div class="annotation">{annotation}</div>

    <div class="conf-row">
      <div class="conf-label">Confidence</div>
      {confidence_bar(tier)}
      <div class="conf-tier" style="color:{confidence_color(tier)}">{tier}</div>
    </div>

    <div class="divider"></div>

    <div>
      <div class="section-label">Evidence chain</div>
      <table>
        <thead>
          <tr>
            <td class="tool-name" style="color:#2a4a63;padding-bottom:4px">tool</td>
            <td class="tool-ann" style="color:#2a4a63;padding-bottom:4px">result</td>
            <td class="tool-tier" style="color:#2a4a63;padding-bottom:4px;font-family:'IBM Plex Mono',monospace;font-size:10px">tier</td>
            <td class="score-val" style="color:#2a4a63;padding-bottom:4px">e-val</td>
          </tr>
        </thead>
        <tbody>
          {evidence_rows}
          {skipped_rows}
        </tbody>
      </table>
      <div class="conf-legend">
        <span><span class="badge" style="background:#4ade8020;border:1px solid #4ade8040;color:#4ade80">HIGH</span> strong homology, named protein</span>
        <span><span class="badge" style="background:#facc1520;border:1px solid #facc1540;color:#facc15">MOD</span> meaningful signal, incomplete</span>
        <span><span class="badge" style="background:#fb923c20;border:1px solid #fb923c40;color:#fb923c">LOW</span> weak or speculative evidence</span>
        <span><span class="badge" style="background:#94a3b820;border:1px solid #94a3b840;color:#94a3b8">UNK</span> no hits above threshold</span>
      </div>
    </div>

    <div class="divider"></div>

    <div class="plddt-section">
      <div class="section-label">Structure confidence (pLDDT)</div>
      <div style="display:flex;align-items:center;gap:0.75rem">
        <div class="plddt-mean">Mean {plddt["mean"]} <span style="font-size:9px;color:#3d5a73;font-weight:400">/ 100</span></div>
        <div style="font-size:10px;color:#3d5a73;font-family:'IBM Plex Mono',monospace">
          {plddt["total"]} residues
        </div>
      </div>
      <div class="plddt-bar-wrap">
        <div style="height:100%;display:flex">{plddt_bar_html}</div>
      </div>
      <div class="plddt-legend">
        <span style="--c:#0053d6"><span style="background:#0053d6;width:8px;height:8px;border-radius:2px;display:inline-block"></span>≥90 very high</span>
        <span><span style="background:#65cbf3;width:8px;height:8px;border-radius:2px;display:inline-block"></span>70–90 high</span>
        <span><span style="background:#ffdb13;width:8px;height:8px;border-radius:2px;display:inline-block"></span>50–70 medium</span>
        <span><span style="background:#ff7d45;width:8px;height:8px;border-radius:2px;display:inline-block"></span>&lt;50 low</span>
      </div>
    </div>

    {warning_html}

    <div class="meta-row" style="margin-top:auto">
      <div class="meta-item"><strong>Length</strong> {seq_len} aa</div>
      <div class="meta-item"><strong>Structure</strong> ESMFold v1</div>
      <div class="meta-item"><strong>Generated</strong> {report.get("generated","")[:10]}</div>
    </div>
  </div>

  <!-- RIGHT: 3D canvas -->
  <div class="right">
    <canvas id="mol"></canvas>
    <div class="structure-label">ESMFold predicted structure · Cα backbone</div>
    <div class="plddt-legend-3d">
      <span style="color:#0053d6"><span style="background:#0053d6"></span>very high</span>
      <span style="color:#65cbf3"><span style="background:#65cbf3"></span>high</span>
      <span style="color:#ffdb13"><span style="background:#ffdb13"></span>medium</span>
      <span style="color:#ff7d45"><span style="background:#ff7d45"></span>low</span>
    </div>
  </div>

</div>

<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
<script>
const ATOMS = {atoms_json};

// pLDDT -> hex color
function plddt2color(p) {{
  if (p >= 90) return 0x0053d6;
  if (p >= 70) return 0x65cbf3;
  if (p >= 50) return 0xffdb13;
  return 0xff7d45;
}}

const canvas = document.getElementById('mol');
const W = canvas.parentElement.clientWidth;
const H = canvas.parentElement.clientHeight;

const renderer = new THREE.WebGLRenderer({{ canvas, antialias: true, alpha: true }});
renderer.setSize(W, H);
renderer.setPixelRatio(window.devicePixelRatio);
renderer.setClearColor(0x000000, 0);

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(45, W / H, 0.1, 1000);
camera.position.set(0, 0, 120);

// Ambient + directional light
scene.add(new THREE.AmbientLight(0xffffff, 0.4));
const dir = new THREE.DirectionalLight(0xffffff, 0.8);
dir.position.set(1, 2, 3);
scene.add(dir);

// Center atoms
const xs = ATOMS.map(a => a[0]);
const ys = ATOMS.map(a => a[1]);
const zs = ATOMS.map(a => a[2]);
const cx = (Math.min(...xs) + Math.max(...xs)) / 2;
const cy = (Math.min(...ys) + Math.max(...ys)) / 2;
const cz = (Math.min(...zs) + Math.max(...zs)) / 2;
const span = Math.max(
  Math.max(...xs) - Math.min(...xs),
  Math.max(...ys) - Math.min(...ys),
  Math.max(...zs) - Math.min(...zs)
);
camera.position.z = span * 1.5;

const group = new THREE.Group();
scene.add(group);

// Draw Cα spheres colored by pLDDT
const sphGeo = new THREE.SphereGeometry(0.6, 8, 8);
ATOMS.forEach(([x, y, z, p]) => {{
  const mat = new THREE.MeshPhongMaterial({{ color: plddt2color(p), shininess: 60 }});
  const sphere = new THREE.Mesh(sphGeo, mat);
  sphere.position.set(x - cx, y - cy, z - cz);
  group.add(sphere);
}});

// Draw backbone bonds between consecutive Cα
for (let i = 0; i < ATOMS.length - 1; i++) {{
  const [x1,y1,z1,p1] = ATOMS[i];
  const [x2,y2,z2,p2] = ATOMS[i+1];
  const dist = Math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2);
  if (dist > 10) continue; // skip chain breaks

  const dir3 = new THREE.Vector3(x2-x1, y2-y1, z2-z1).normalize();
  const mid = new THREE.Vector3((x1+x2)/2-cx, (y1+y2)/2-cy, (z1+z2)/2-cz);

  const cyl = new THREE.CylinderGeometry(0.25, 0.25, dist, 6);
  const mat = new THREE.MeshPhongMaterial({{ color: plddt2color((p1+p2)/2), shininess: 40 }});
  const mesh = new THREE.Mesh(cyl, mat);
  mesh.position.copy(mid);
  mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0,1,0), dir3);
  group.add(mesh);
}}

// Animate: slow spin
let angle = 0;
function animate() {{
  requestAnimationFrame(animate);
  angle += 0.008;
  group.rotation.y = angle;
  group.rotation.x = Math.sin(angle * 0.3) * 0.15;
  renderer.render(scene, camera);
}}
animate();

// Resize
window.addEventListener('resize', () => {{
  const w = canvas.parentElement.clientWidth;
  const h = canvas.parentElement.clientHeight;
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
  renderer.setSize(w, h);
}});
</script>
</body>
</html>"""


def main():
    parser = argparse.ArgumentParser(
        prog="generate_card",
        description="Generate an animated protein annotation card from an ORACLE report."
    )
    parser.add_argument("report", help="Path to ORACLE .report.json file")
    parser.add_argument(
        "--fasta", metavar="FASTA",
        help="FASTA file containing the sequence (alternative to --sequence)"
    )
    parser.add_argument(
        "--sequence", metavar="SEQ",
        help="Raw amino acid sequence string"
    )
    parser.add_argument(
        "-o", "--output", metavar="HTML",
        help="Output HTML file path. Default: same name as report with .card.html"
    )
    args = parser.parse_args()

    # Load report
    print(f"Loading report: {args.report}")
    report = load_report(args.report)

    # Get sequence
    if args.fasta:
        _, sequence = load_sequence(args.fasta)
    elif args.sequence:
        sequence = args.sequence.strip().replace("\n", "").replace(" ", "")
    else:
        print("Error: provide --fasta or --sequence")
        sys.exit(1)

    if len(sequence) > 400:
        print(f"Sequence is {len(sequence)} aa, which exceeds the ESMFold API limit (~400 aa).")
        print()
        print("The card generator is intended for uncharacterized or dark matter sequences")
        print("under 400 aa where ESMFold structure prediction is the only available")
        print("structural hypothesis. For longer well-characterized proteins, the")
        print("annotation report is the primary output.")
        print()
        print("If you need a card for this sequence, consider:")
        print("  - Using the annotation report directly (reports/*.report.txt)")
        print("  - Fetching the experimental PDB structure if one exists")
        sys.exit(0)

    # Predict structure
    print("Calling ESMFold API...")
    try:
        pdb_str = fold_sequence(sequence)
        plddt = compute_plddt_stats(pdb_str)
        print(f"  Structure predicted. Mean pLDDT: {plddt['mean']} ({plddt['total']} residues)")
    except Exception as e:
        print(f"ESMFold failed: {e}")
        print("Generating card without structure...")
        pdb_str = ""
        plddt = {"mean": 0, "very_high": 0, "high": 0, "medium": 0, "low": 0, "total": 0}

    # Determine output path
    if args.output:
        out_path = args.output
    else:
        base = args.report.replace(".report.json", "").replace(".json", "")
        out_path = f"{base}.card.html"

    # Generate card
    print("Generating card...")
    html = generate_html(report, pdb_str, plddt)

    os.makedirs(os.path.dirname(out_path) if os.path.dirname(out_path) else ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"Card saved: {out_path}")
    print(f"Open in a browser to view the animated structure.")


if __name__ == "__main__":
    main()