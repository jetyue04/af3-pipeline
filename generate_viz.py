"""Generate GTP and GDP complex HTML visualizations with interactive distance filters."""

import json
from pathlib import Path

RESULTS = Path("/Users/jetyue04/af3/af3-pipeline/results")
RESULTS.mkdir(parents=True, exist_ok=True)

CONFIGS = {
    "gtp": {
        "cif":    "/Users/jetyue04/af3/output/fam129b_rab5a_gtp_mg_complex/fam129b_rab5a_gtp_mg_complex_model.cif",
        "scores": "/Users/jetyue04/af3/output/fam129b_rab5a_gtp_mg_complex/fam129b_rab5a_gtp_mg_complex_summary_confidences.json",
        "ligand_label": "GTP",
        "ligand_chain": "C",
        "contacts": {
            4: {
                "rab5a": [(22,'K','P-loop'),(53,'I','other'),(54,'G','other'),(55,'A','other'),
                          (56,'A','other'),(57,'F','SwitchII'),(59,'T','SwitchII'),(70,'K','SwitchII'),
                          (74,'W','SwitchII'),(81,'R','other'),(82,'Y','other'),(85,'L','other'),
                          (88,'M','other'),(89,'Y','other')],
                "fam":   [(359,'T'),(362,'R'),(363,'D'),(366,'F'),(367,'K'),(370,'T'),(371,'D'),
                          (373,'N'),(374,'L'),(377,'I'),(463,'E'),(466,'C'),(467,'K'),(470,'Q')],
            },
            5: {
                "rab5a": [(20,'Q','P-loop'),(22,'K','P-loop'),(50,'E','other'),(53,'I','other'),
                          (54,'G','other'),(55,'A','other'),(56,'A','other'),(57,'F','SwitchII'),
                          (58,'L','SwitchII'),(59,'T','SwitchII'),(70,'K','SwitchII'),(72,'E','SwitchII'),
                          (74,'W','SwitchII'),(81,'R','other'),(82,'Y','other'),(85,'L','other'),
                          (88,'M','other'),(89,'Y','other')],
                "fam":   [(359,'T'),(362,'R'),(363,'D'),(366,'F'),(367,'K'),(370,'T'),(371,'D'),
                          (373,'N'),(374,'L'),(377,'I'),(378,'N'),(463,'E'),(466,'C'),(467,'K'),
                          (470,'Q'),(474,'E'),(477,'L'),(478,'K')],
            },
            8: {
                "rab5a": [(20,'Q','P-loop'),(22,'K','P-loop'),(41,'V','SwitchI'),(42,'K','SwitchI'),
                          (50,'E','other'),(51,'S','other'),(52,'T','other'),(53,'I','other'),
                          (54,'G','other'),(55,'A','other'),(56,'A','other'),(57,'F','SwitchII'),
                          (58,'L','SwitchII'),(59,'T','SwitchII'),(60,'Q','SwitchII'),(61,'T','SwitchII'),
                          (70,'K','SwitchII'),(71,'F','SwitchII'),(72,'E','SwitchII'),(73,'I','SwitchII'),
                          (74,'W','SwitchII'),(75,'D','other'),(76,'T','other'),(77,'A','other'),
                          (81,'R','other'),(82,'Y','other'),(85,'L','other'),(88,'M','other'),(89,'Y','other')],
                "fam":   [(356,'Q'),(359,'T'),(360,'E'),(362,'R'),(363,'D'),(364,'V'),(366,'F'),
                          (367,'K'),(368,'E'),(369,'V'),(370,'T'),(371,'D'),(372,'M'),(373,'N'),
                          (374,'L'),(375,'N'),(377,'I'),(378,'N'),(461,'T'),(462,'K'),(463,'E'),
                          (464,'E'),(465,'L'),(466,'C'),(467,'K'),(468,'S'),(469,'I'),(470,'Q'),
                          (471,'R'),(473,'L'),(474,'E'),(477,'L'),(478,'K'),(481,'D')],
            },
        },
    },
    "gdp": {
        "cif":    "/Users/jetyue04/af3/output/fam129b_rab5a_gdp_mg_complex/fam129b_rab5a_gdp_mg_complex_model.cif",
        "scores": "/Users/jetyue04/af3/output/fam129b_rab5a_gdp_mg_complex/fam129b_rab5a_gdp_mg_complex_summary_confidences.json",
        "ligand_label": "GDP",
        "ligand_chain": "C",
        "contacts": {
            4: {
                "rab5a": [(22,'K','P-loop'),(50,'E','other'),(53,'I','other'),(54,'G','other'),
                          (55,'A','other'),(56,'A','other'),(57,'F','SwitchII'),(59,'T','SwitchII'),
                          (70,'K','SwitchII'),(72,'E','SwitchII'),(74,'W','SwitchII'),(81,'R','other'),
                          (82,'Y','other'),(85,'L','other'),(88,'M','other'),(89,'Y','other')],
                "fam":   [(362,'R'),(363,'D'),(366,'F'),(367,'K'),(370,'T'),(371,'D'),(373,'N'),
                          (374,'L'),(377,'I'),(463,'E'),(466,'C'),(467,'K'),(470,'Q'),(474,'E'),
                          (477,'L'),(478,'K')],
            },
            5: {
                "rab5a": [(22,'K','P-loop'),(50,'E','other'),(52,'T','other'),(53,'I','other'),
                          (54,'G','other'),(55,'A','other'),(56,'A','other'),(57,'F','SwitchII'),
                          (58,'L','SwitchII'),(59,'T','SwitchII'),(70,'K','SwitchII'),(72,'E','SwitchII'),
                          (74,'W','SwitchII'),(81,'R','other'),(82,'Y','other'),(85,'L','other'),
                          (88,'M','other'),(89,'Y','other')],
                "fam":   [(359,'T'),(362,'R'),(363,'D'),(366,'F'),(367,'K'),(370,'T'),(371,'D'),
                          (373,'N'),(374,'L'),(377,'I'),(463,'E'),(466,'C'),(467,'K'),(470,'Q'),
                          (474,'E'),(477,'L'),(478,'K'),(481,'D')],
            },
            8: {
                "rab5a": [(20,'Q','P-loop'),(22,'K','P-loop'),(34,'S','SwitchI'),(41,'V','SwitchI'),
                          (42,'K','SwitchI'),(50,'E','other'),(51,'S','other'),(52,'T','other'),
                          (53,'I','other'),(54,'G','other'),(55,'A','other'),(56,'A','other'),
                          (57,'F','SwitchII'),(58,'L','SwitchII'),(59,'T','SwitchII'),(60,'Q','SwitchII'),
                          (61,'T','SwitchII'),(70,'K','SwitchII'),(71,'F','SwitchII'),(72,'E','SwitchII'),
                          (73,'I','SwitchII'),(74,'W','SwitchII'),(75,'D','other'),(76,'T','other'),
                          (81,'R','other'),(82,'Y','other'),(85,'L','other'),(88,'M','other'),(89,'Y','other')],
                "fam":   [(359,'T'),(360,'E'),(362,'R'),(363,'D'),(364,'V'),(365,'F'),(366,'F'),
                          (367,'K'),(368,'E'),(369,'V'),(370,'T'),(371,'D'),(372,'M'),(373,'N'),
                          (374,'L'),(375,'N'),(376,'V'),(377,'I'),(378,'N'),(461,'T'),(462,'K'),
                          (463,'E'),(464,'E'),(465,'L'),(466,'C'),(467,'K'),(468,'S'),(469,'I'),
                          (470,'Q'),(471,'R'),(473,'L'),(474,'E'),(475,'R'),(476,'V'),(477,'L'),
                          (478,'K'),(479,'K'),(481,'D')],
            },
        },
    },
}

REGION_COLOR = {"P-loop": "#facc15", "SwitchI": "#fb923c", "SwitchII": "#f87171", "other": "#a78bfa"}
REGION_LABEL = {"P-loop": "P-loop (grips phosphate)", "SwitchI": "Switch I", "SwitchII": "Switch II (state sensor)", "other": "Other"}


def resi_str(lst):
    return ",".join(str(x[0]) for x in lst)


def make_table_rows(contacts, side):
    rows = []
    for item in contacts:
        num, aa = item[0], item[1]
        reg = item[2] if len(item) > 2 else "other"
        color = REGION_COLOR[reg]
        label = REGION_LABEL[reg]
        rows.append(f'<tr><td>{num}</td><td style="font-weight:600">{aa}</td>'
                    f'<td><span class="tag" style="background:{color}22;color:{color};border:1px solid {color}55">{label}</span></td></tr>')
    return "\n".join(rows)


def build_html(state, cfg):
    cif_text = Path(cfg["cif"]).read_text()
    cif_escaped = cif_text.replace("\\", "\\\\").replace("`", "\\`")

    with open(cfg["scores"]) as f:
        sc = json.load(f)

    iptm       = sc["iptm"]
    ptm        = sc["ptm"]
    pair_iptm  = sc["chain_pair_iptm"][0][1]
    ranking    = sc["ranking_score"]
    frac_dis   = sc["fraction_disordered"]

    def badge(v, hi=0.75, lo=0.5):
        cls = "good" if v >= hi else "mid" if v >= lo else "low"
        return f'<span class="val {cls}">{v:.2f}</span>'

    contacts = cfg["contacts"]
    lig = cfg["ligand_label"]
    lig_chain = cfg["ligand_chain"]

    # precompute JS data for all cutoffs
    js_data = {}
    for cut, data in contacts.items():
        js_data[cut] = {
            "rab5a_resi": resi_str(data["rab5a"]),
            "fam_resi":   resi_str(data["fam"]),
            "rab5a_rows": make_table_rows(data["rab5a"], "rab5a"),
            "fam_rows":   make_table_rows(data["fam"], "fam"),
            "rab5a_count": len(data["rab5a"]),
            "fam_count":   len(data["fam"]),
        }

    title = f"FAM129B – RAB5A ({lig}-bound) · AlphaFold 3"
    state_color = "#4ade80" if state == "gtp" else "#94a3b8"
    lig_color   = "#fb923c" if state == "gtp" else "#60a5fa"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>{title}</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
* {{ box-sizing:border-box; margin:0; padding:0; }}
body {{ font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; background:#0f1117; color:#e0e0e0; }}
header {{ padding:16px 28px 12px; border-bottom:1px solid #2a2a3a; display:flex; align-items:center; gap:14px; }}
h1 {{ font-size:1.25rem; font-weight:600; color:#fff; }}
.state-badge {{ padding:3px 10px; border-radius:20px; font-size:0.78rem; font-weight:700;
                background:{state_color}22; color:{state_color}; border:1px solid {state_color}55; }}
.layout {{ display:flex; height:calc(100vh - 58px); }}
#viewer {{ flex:1; }}
.panel {{ width:320px; background:#16181f; border-left:1px solid #2a2a3a; overflow-y:auto;
          padding:16px; display:flex; flex-direction:column; gap:14px; }}
.card {{ background:#1e2030; border-radius:8px; padding:14px; }}
.card h3 {{ font-size:0.7rem; text-transform:uppercase; letter-spacing:0.08em; color:#555; margin-bottom:10px; }}
.metric {{ display:flex; justify-content:space-between; align-items:center; margin-bottom:6px; font-size:0.85rem; }}
.val {{ font-weight:700; font-size:0.95rem; }}
.good {{ color:#4ade80; }} .mid {{ color:#facc15; }} .low {{ color:#f87171; }}
.btn-row {{ display:flex; gap:6px; flex-wrap:wrap; }}
button {{ padding:6px 12px; border:none; border-radius:6px; background:#2a2d3e; color:#ccc;
          cursor:pointer; font-size:0.8rem; transition:background 0.15s; }}
button:hover {{ background:#363a50; }}
button.active {{ background:#3b4fd8; color:#fff; }}
.legend-item {{ display:flex; align-items:center; gap:8px; margin-bottom:7px; font-size:0.83rem; }}
.swatch {{ width:12px; height:12px; border-radius:3px; flex-shrink:0; }}
table {{ width:100%; border-collapse:collapse; font-size:0.78rem; }}
th {{ text-align:left; padding:4px 6px; color:#555; font-weight:600; border-bottom:1px solid #2a2a3a; }}
td {{ padding:4px 6px; border-bottom:1px solid #1a1c28; }}
.tag {{ padding:1px 6px; border-radius:4px; font-size:0.7rem; white-space:nowrap; }}
.note {{ font-size:0.78rem; color:#777; line-height:1.55; }}
.note b {{ color:#aaa; }}
.warn {{ color:#f87171; font-size:0.78rem; margin-top:6px; line-height:1.5; }}
</style>
</head>
<body>
<header>
  <h1>FAM129B – RAB5A Complex</h1>
  <div class="state-badge">{lig}-bound</div>
</header>
<div class="layout">
  <div id="viewer"></div>
  <div class="panel">

    <div class="card">
      <h3>AF3 Confidence Scores</h3>
      <div class="metric"><span>Overall iptm</span>{badge(iptm)}</div>
      <div class="metric"><span>ptm</span>{badge(ptm, 0.7, 0.5)}</div>
      <div class="metric"><span>FAM129B↔RAB5A iptm</span>{badge(pair_iptm)}</div>
      <div class="metric"><span>Ranking score</span>{badge(ranking, 0.6, 0.4)}</div>
      <div class="metric"><span>Disordered</span><span class="val {'good' if frac_dis<0.2 else 'mid' if frac_dis<0.35 else 'low'}">{frac_dis:.0%}</span></div>
      <p class="warn">⚠ FAM129B↔RAB5A pair iptm = {pair_iptm:.2f} — AF3 is <b>not confident</b> in this interface.
      {'GDP scores higher than GTP — opposite of what a true GTP-state effector should show.' if state=='gtp' else 'This scores higher than GTP, which is unexpected for a GTP-state effector.'}</p>
    </div>

    <div class="card">
      <h3>Distance Filter — contacts shown</h3>
      <div class="btn-row">
        <button id="cut4" onclick="setCutoff(4)">4 Å — direct</button>
        <button id="cut5" class="active" onclick="setCutoff(5)">5 Å — standard</button>
        <button id="cut8" onclick="setCutoff(8)">8 Å — neighborhood</button>
      </div>
      <p class="note" style="margin-top:8px">
        <b style="color:#e0e0e0">&lt;4 Å</b> = atoms nearly touching (H-bond / vdW range)<br>
        <b style="color:#e0e0e0">4–5 Å</b> = close contact, likely interacting<br>
        <b style="color:#e0e0e0">5–8 Å</b> = nearby but not necessarily bonded
      </p>
    </div>

    <div class="card">
      <h3>View Mode</h3>
      <div class="btn-row">
        <button onclick="showCartoon()">Cartoon</button>
        <button onclick="showSurface()">Surface</button>
        <button onclick="showPLDDT()">pLDDT</button>
        <button onclick="zoomInterface()">Zoom interface</button>
        <button onclick="resetView()">Reset</button>
      </div>
    </div>

    <div class="card">
      <h3>Chain Legend</h3>
      <div class="legend-item"><div class="swatch" style="background:#60a5fa"></div>FAM129B (chain A)</div>
      <div class="legend-item"><div class="swatch" style="background:#4ade80"></div>RAB5A (chain B)</div>
      <div class="legend-item"><div class="swatch" style="background:{lig_color}"></div>{lig} (chain {lig_chain})</div>
      <div class="legend-item"><div class="swatch" style="background:#e879f9"></div>Mg²⁺ (chain D)</div>
      <hr style="border-color:#2a2a3a;margin:8px 0">
      <div class="legend-item"><div class="swatch" style="background:#facc15"></div>P-loop (grips phosphate)</div>
      <div class="legend-item"><div class="swatch" style="background:#fb923c"></div>Switch I</div>
      <div class="legend-item"><div class="swatch" style="background:#f87171"></div>Switch II (state sensor)</div>
      <div class="legend-item"><div class="swatch" style="background:#a78bfa"></div>Other interface residues</div>
    </div>

    <div class="card">
      <h3>RAB5A interface residues — <span id="rab5a-count">{js_data[5]['rab5a_count']}</span> at current cutoff</h3>
      <table>
        <tr><th>#</th><th>AA</th><th>Region</th></tr>
        <tbody id="rab5a-table">{js_data[5]['rab5a_rows']}</tbody>
      </table>
    </div>

    <div class="card">
      <h3>FAM129B interface residues — <span id="fam-count">{js_data[5]['fam_count']}</span> at current cutoff</h3>
      <table>
        <tr><th>#</th><th>AA</th><th>Region</th></tr>
        <tbody id="fam-table">{js_data[5]['fam_rows']}</tbody>
      </table>
    </div>

    <div class="card">
      <h3>What low scores mean</h3>
      <p class="note">
        <b>iptm</b> measures how confident AF3 is in the <i>relative position</i> of two chains.
        Scores &gt;0.75 = strong predicted binding. Scores &lt;0.4 = uncertain.
        <br><br>
        This complex scores {iptm:.2f}. Possible reasons:
        <br><br>
        <b>1.</b> FAM129B may bind RAB5A <i>indirectly</i> (through another protein)<br>
        <b>2.</b> The interaction may be transient or weak — AF3 struggles with these<br>
        <b>3.</b> FAM129B may need a post-translational modification to bind<br>
        <b>4.</b> Only a specific domain of FAM129B (not what we modelled) may be responsible
        <br><br>
        The contacts landing on <b>Switch II</b> (the GTP state sensor) are still biologically
        meaningful — just treat the pose geometry with caution until validated experimentally.
      </p>
    </div>

  </div>
</div>

<script>
const CIF = `{cif_escaped}`;

const CONTACTS = {{
  4:  {{ rab5a: "{js_data[4]['rab5a_resi']}", fam: "{js_data[4]['fam_resi']}",
         rab5aRows: `{js_data[4]['rab5a_rows']}`, famRows: `{js_data[4]['fam_rows']}`,
         rab5aCount: {js_data[4]['rab5a_count']}, famCount: {js_data[4]['fam_count']} }},
  5:  {{ rab5a: "{js_data[5]['rab5a_resi']}", fam: "{js_data[5]['fam_resi']}",
         rab5aRows: `{js_data[5]['rab5a_rows']}`, famRows: `{js_data[5]['fam_rows']}`,
         rab5aCount: {js_data[5]['rab5a_count']}, famCount: {js_data[5]['fam_count']} }},
  8:  {{ rab5a: "{js_data[8]['rab5a_resi']}", fam: "{js_data[8]['fam_resi']}",
         rab5aRows: `{js_data[8]['rab5a_rows']}`, famRows: `{js_data[8]['fam_rows']}`,
         rab5aCount: {js_data[8]['rab5a_count']}, famCount: {js_data[8]['fam_count']} }},
}};

const REGION_COLORS = {{
  "P-loop":   "#facc15",
  "SwitchI":  "#fb923c",
  "SwitchII": "#f87171",
  "other":    "#a78bfa",
}};

// RAB5A residues with their regions for colouring
const RAB5A_REGIONS = {{
  "P-loop":   [17,18,19,20,21,22,23,24],
  "SwitchI":  [33,34,35,36,37,38,39,40,41,42],
  "SwitchII": [57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74],
}};

let currentCutoff = 5;
let currentView = "cartoon";

const viewer = $3Dmol.createViewer("viewer", {{ backgroundColor: "0x0f1117" }});
viewer.addModel(CIF, "cif");

function applyBase() {{
  viewer.setStyle({{}}, {{}});
  viewer.setStyle({{chain:"A"}}, {{cartoon:{{color:"#60a5fa"}}}});
  viewer.setStyle({{chain:"B"}}, {{cartoon:{{color:"#4ade80"}}}});
  viewer.setStyle({{chain:"{lig_chain}"}}, {{stick:{{colorscheme:"{('orange' if state=='gtp' else 'lightblue')}Carbon", radius:0.28}}}});
  viewer.setStyle({{chain:"D"}}, {{sphere:{{color:"#e879f9", radius:0.65}}}});
}}

function applyRegionColors() {{
  for (const [region, residues] of Object.entries(RAB5A_REGIONS)) {{
    viewer.setStyle({{chain:"B", resi:residues.join(",")}},
      {{cartoon:{{color:REGION_COLORS[region]}}, stick:{{colorscheme:"greenCarbon", radius:0.0}}}});
  }}
}}

function applyContacts(cutoff) {{
  const d = CONTACTS[cutoff];
  // FAM129B contact residues — cyan sticks
  viewer.setStyle({{chain:"A", resi:d.fam}},   {{cartoon:{{color:"#60a5fa"}}, stick:{{color:"#22d3ee", radius:0.22}}}});
  // RAB5A contact residues — coloured by region
  const rab5aList = d.rab5a.split(",").map(Number);
  for (const resnum of rab5aList) {{
    let col = REGION_COLORS["other"];
    for (const [region, residues] of Object.entries(RAB5A_REGIONS)) {{
      if (residues.includes(resnum)) {{ col = REGION_COLORS[region]; break; }}
    }}
    viewer.setStyle({{chain:"B", resi:String(resnum)}}, {{cartoon:{{color:col}}, stick:{{color:col, radius:0.22}}}});
  }}
  // Labels on RAB5A contacts
  viewer.removeAllLabels();
  for (const resnum of rab5aList) {{
    viewer.addResLabels({{chain:"B", resi:String(resnum)}},
      {{fontSize:10, fontColor:"white", showBackground:false, backgroundOpacity:0}});
  }}
}}

function updateTable(cutoff) {{
  const d = CONTACTS[cutoff];
  document.getElementById("rab5a-table").innerHTML = d.rab5aRows;
  document.getElementById("fam-table").innerHTML   = d.famRows;
  document.getElementById("rab5a-count").textContent = d.rab5aCount;
  document.getElementById("fam-count").textContent   = d.famCount;
  [4,5,8].forEach(c => document.getElementById("cut"+c).classList.remove("active"));
  document.getElementById("cut"+cutoff).classList.add("active");
}}

function setCutoff(cutoff) {{
  currentCutoff = cutoff;
  applyBase();
  if (currentView === "surface") {{
    viewer.addSurface($3Dmol.SurfaceType.VDW, {{color:"#60a5fa", opacity:0.4}}, {{chain:"A"}});
    viewer.setStyle({{chain:"A"}}, {{cartoon:{{color:"#60a5fa", opacity:0.2}}}});
  }} else if (currentView === "plddt") {{
    viewer.setStyle({{chain:"A"}}, {{cartoon:{{colorfunc: plddt}}}});
    viewer.setStyle({{chain:"B"}}, {{cartoon:{{colorfunc: plddt}}}});
  }} else {{
    applyRegionColors();
  }}
  applyContacts(cutoff);
  updateTable(cutoff);
  viewer.render();
}}

function showCartoon() {{
  currentView = "cartoon";
  viewer.removeAllSurfaces();
  applyBase();
  applyRegionColors();
  applyContacts(currentCutoff);
  viewer.render();
}}

function showSurface() {{
  currentView = "surface";
  viewer.removeAllSurfaces();
  applyBase();
  viewer.setStyle({{chain:"A"}}, {{cartoon:{{color:"#60a5fa", opacity:0.2}}}});
  viewer.addSurface($3Dmol.SurfaceType.VDW, {{color:"#60a5fa", opacity:0.4}}, {{chain:"A"}});
  applyRegionColors();
  applyContacts(currentCutoff);
  viewer.render();
}}

function showPLDDT() {{
  currentView = "plddt";
  viewer.removeAllSurfaces();
  viewer.setStyle({{}}, {{}});
  viewer.setStyle({{chain:"A"}}, {{cartoon:{{colorfunc: plddt}}}});
  viewer.setStyle({{chain:"B"}}, {{cartoon:{{colorfunc: plddt}}}});
  viewer.setStyle({{chain:"{lig_chain}"}}, {{stick:{{colorscheme:"orangeCarbon", radius:0.28}}}});
  viewer.setStyle({{chain:"D"}}, {{sphere:{{color:"#e879f9", radius:0.65}}}});
  applyContacts(currentCutoff);
  viewer.render();
}}

function zoomInterface() {{
  const d = CONTACTS[currentCutoff];
  viewer.zoomTo({{chain:"B", resi:d.rab5a}});
  viewer.render();
}}

function resetView() {{
  viewer.zoomTo();
  showCartoon();
}}

function plddt(atom) {{
  const b = atom.b;
  if (b >= 90) return "#2563eb";
  if (b >= 70) return "#60a5fa";
  if (b >= 50) return "#facc15";
  return "#f87171";
}}

showCartoon();
viewer.zoomTo();
viewer.render();
</script>
</body>
</html>"""

    return html


for state, cfg in CONFIGS.items():
    html = build_html(state, cfg)
    out = RESULTS / f"fam129b_rab5a_{state}_complex.html"
    out.write_text(html)
    print(f"Saved → {out}")
