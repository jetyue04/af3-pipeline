"""Generate a shareable HTML visualization of the FAM129B–RAB5A–GTP complex."""

import json
from pathlib import Path

CIF = Path("/Users/jetyue04/af3/output/fam129b_rab5a_gtp_mg_complex/fam129b_rab5a_gtp_mg_complex_model.cif")
SCORES = Path("/Users/jetyue04/af3/output/fam129b_rab5a_gtp_mg_complex/fam129b_rab5a_gtp_mg_complex_summary_confidences.json")
OUT = Path("/Users/jetyue04/af3/af3-pipeline/results/fam129b_rab5a_gtp_complex.html")

OUT.parent.mkdir(parents=True, exist_ok=True)

cif_text = CIF.read_text()

with open(SCORES) as f:
    scores = json.load(f)

iptm = scores["iptm"]
ptm  = scores["ptm"]
ranking = scores["ranking_score"]
frac_dis = scores["fraction_disordered"]

# chain_pair_iptm[0][1] = FAM129B vs RAB5A confidence
pair_iptm = scores["chain_pair_iptm"][0][1]

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>FAM129B – RAB5A – GTP Complex (AlphaFold 3)</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         background: #0f1117; color: #e0e0e0; }}
  header {{ padding: 20px 32px 12px; border-bottom: 1px solid #2a2a3a; }}
  h1 {{ font-size: 1.4rem; font-weight: 600; color: #fff; }}
  .subtitle {{ font-size: 0.85rem; color: #888; margin-top: 4px; }}
  .layout {{ display: flex; height: calc(100vh - 80px); }}
  #viewer {{ flex: 1; }}
  .panel {{ width: 280px; background: #16181f; border-left: 1px solid #2a2a3a;
            padding: 20px; overflow-y: auto; display: flex; flex-direction: column; gap: 18px; }}
  .card {{ background: #1e2030; border-radius: 8px; padding: 14px; }}
  .card h3 {{ font-size: 0.75rem; text-transform: uppercase; letter-spacing: 0.08em;
              color: #666; margin-bottom: 10px; }}
  .metric {{ display: flex; justify-content: space-between; align-items: center;
             margin-bottom: 6px; font-size: 0.88rem; }}
  .metric .val {{ font-weight: 600; font-size: 1rem; }}
  .good  {{ color: #4ade80; }}
  .mid   {{ color: #facc15; }}
  .low   {{ color: #f87171; }}
  .legend-item {{ display: flex; align-items: center; gap: 8px;
                  margin-bottom: 8px; font-size: 0.88rem; }}
  .swatch {{ width: 14px; height: 14px; border-radius: 3px; flex-shrink: 0; }}
  .controls button {{
    width: 100%; padding: 8px; margin-bottom: 6px; border: none; border-radius: 6px;
    background: #2a2d3e; color: #ccc; cursor: pointer; font-size: 0.85rem;
    transition: background 0.15s;
  }}
  .controls button:hover {{ background: #363a50; }}
  .note {{ font-size: 0.78rem; color: #666; line-height: 1.5; }}
</style>
</head>
<body>

<header>
  <h1>FAM129B – RAB5A (GTP-bound) Complex</h1>
  <div class="subtitle">AlphaFold 3 prediction · seed-1 sample-0 (top-ranked model)</div>
</header>

<div class="layout">
  <div id="viewer"></div>

  <div class="panel">

    <div class="card">
      <h3>Confidence Scores</h3>
      <div class="metric">
        <span>Overall iptm</span>
        <span class="val {'good' if iptm >= 0.75 else 'mid' if iptm >= 0.5 else 'low'}">{iptm:.2f}</span>
      </div>
      <div class="metric">
        <span>ptm</span>
        <span class="val {'good' if ptm >= 0.7 else 'mid' if ptm >= 0.5 else 'low'}">{ptm:.2f}</span>
      </div>
      <div class="metric">
        <span>FAM129B↔RAB5A iptm</span>
        <span class="val {'good' if pair_iptm >= 0.75 else 'mid' if pair_iptm >= 0.5 else 'low'}">{pair_iptm:.2f}</span>
      </div>
      <div class="metric">
        <span>Ranking score</span>
        <span class="val mid">{ranking:.2f}</span>
      </div>
      <div class="metric">
        <span>Disordered</span>
        <span class="val {'good' if frac_dis < 0.2 else 'mid' if frac_dis < 0.35 else 'low'}">{frac_dis:.0%}</span>
      </div>
    </div>

    <div class="card">
      <h3>Chain Legend</h3>
      <div class="legend-item"><div class="swatch" style="background:#60a5fa"></div>FAM129B (chain A)</div>
      <div class="legend-item"><div class="swatch" style="background:#4ade80"></div>RAB5A (chain B)</div>
      <div class="legend-item"><div class="swatch" style="background:#fb923c"></div>GTP (chain C)</div>
      <div class="legend-item"><div class="swatch" style="background:#e879f9"></div>Mg²⁺ (chain D)</div>
    </div>

    <div class="card controls">
      <h3>Views</h3>
      <button onclick="showCartoon()">Cartoon (overview)</button>
      <button onclick="showSurface()">Surface (FAM129B) + Cartoon (RAB5A)</button>
      <button onclick="showInterface()">Interface close-up</button>
      <button onclick="showPLDDT()">pLDDT confidence colouring</button>
      <button onclick="resetView()">Reset view</button>
    </div>

    <div class="card">
      <h3>What you're looking at</h3>
      <p class="note">
        <b style="color:#4ade80">RAB5A</b> is a molecular switch that controls how cells sort
        cargo inside vesicles. When loaded with GTP (orange) it's "switched on."
        <br><br>
        <b style="color:#60a5fa">FAM129B</b> is a protein we think docks onto RAB5A when it's
        active — this could affect how certain receptors (like PD-L1) are recycled in cancer cells.
        <br><br>
        AlphaFold 3 predicted this 3D shape from the protein sequences alone.
        We're now designing small molecules that could block this interaction.
      </p>
    </div>

  </div>
</div>

<script>
const CIF = `{cif_text.replace(chr(96), chr(92) + chr(96))}`;

const viewer = $3Dmol.createViewer("viewer", {{
  backgroundColor: "0x0f1117"
}});

viewer.addModel(CIF, "cif");

function showCartoon() {{
  viewer.setStyle({{}}, {{}});
  viewer.setStyle({{chain:"A"}}, {{cartoon:{{color:"#60a5fa", opacity:1}}}});
  viewer.setStyle({{chain:"B"}}, {{cartoon:{{color:"#4ade80", opacity:1}}}});
  viewer.setStyle({{chain:"C"}}, {{stick:{{colorscheme:"orangeCarbon", radius:0.25}}}});
  viewer.setStyle({{chain:"D"}}, {{sphere:{{color:"#e879f9", radius:0.6}}}});
  viewer.render();
}}

function showSurface() {{
  viewer.setStyle({{}}, {{}});
  viewer.setStyle({{chain:"A"}}, {{cartoon:{{color:"#60a5fa", opacity:0.3}}}});
  viewer.addSurface($3Dmol.SurfaceType.VDW, {{color:"#60a5fa", opacity:0.45}}, {{chain:"A"}});
  viewer.setStyle({{chain:"B"}}, {{cartoon:{{color:"#4ade80"}}}});
  viewer.setStyle({{chain:"C"}}, {{stick:{{colorscheme:"orangeCarbon", radius:0.3}}}});
  viewer.setStyle({{chain:"D"}}, {{sphere:{{color:"#e879f9", radius:0.6}}}});
  viewer.render();
}}

function showInterface() {{
  viewer.setStyle({{}}, {{}});
  viewer.setStyle({{chain:"A"}}, {{cartoon:{{color:"#60a5fa", opacity:0.25}}}});
  viewer.setStyle({{chain:"B"}}, {{cartoon:{{color:"#4ade80", opacity:0.25}}}});
  // highlight residues within ~8Å of the other chain (approximate Switch I/II region)
  viewer.setStyle({{chain:"A", resi:"1-50"}},  {{cartoon:{{color:"#60a5fa", opacity:0.8}}}});
  viewer.setStyle({{chain:"B", resi:"30-60"}}, {{cartoon:{{color:"#4ade80", opacity:0.8}}, stick:{{colorscheme:"greenCarbon", radius:0.2}}}});
  viewer.setStyle({{chain:"C"}}, {{stick:{{colorscheme:"orangeCarbon", radius:0.35}}}});
  viewer.setStyle({{chain:"D"}}, {{sphere:{{color:"#e879f9", radius:0.7}}}});
  viewer.zoomTo({{chain:"C"}});
  viewer.render();
}}

function showPLDDT() {{
  viewer.setStyle({{}}, {{}});
  viewer.setStyle({{chain:"A"}}, {{cartoon:{{colorfunc: plddt}}}});
  viewer.setStyle({{chain:"B"}}, {{cartoon:{{colorfunc: plddt}}}});
  viewer.setStyle({{chain:"C"}}, {{stick:{{colorscheme:"orangeCarbon", radius:0.25}}}});
  viewer.setStyle({{chain:"D"}}, {{sphere:{{color:"#e879f9", radius:0.6}}}});
  viewer.render();
}}

function plddt(atom) {{
  const b = atom.b; // pLDDT stored in B-factor
  if (b >= 90) return "#2563eb";
  if (b >= 70) return "#60a5fa";
  if (b >= 50) return "#facc15";
  return "#f87171";
}}

function resetView() {{
  showCartoon();
  viewer.zoomTo();
  viewer.render();
}}

showCartoon();
viewer.zoomTo();
viewer.render();
</script>
</body>
</html>
"""

OUT.write_text(html)
print(f"Saved → {OUT}")
print(f"\nScores: iptm={iptm:.2f}  ptm={ptm:.2f}  FAM129B↔RAB5A pair iptm={pair_iptm:.2f}  ranking={ranking:.2f}")
