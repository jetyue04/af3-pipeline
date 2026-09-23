# FAM129B Screening Pipeline — Agreed Plan

**Goal:** Identify small-molecule competitive inhibitors of the RAB5A–FAM129B
protein-protein interface. AF3 co-folding maps the interface; AutoDock Vina +
GNINA screen a PPI-focused compound library against it.

---

## Pipeline Stages

### Stage 1 — AF3 Co-folding (SLURM cluster)
- Input: `protein/inputs/*.json` + SLURM scripts in `jobs/`
- Output: 5 structural samples per job in `AF3_Output/<job_name>/`
- All 11 new jobs ready to submit when on VPN (see `jobs/` for scripts)

**Jobs queued:**

| Job | Type | Notes |
|-----|------|-------|
| `rab7a_monomer` | monomer | |
| `rab7a_full_gdp` | nucleotide state | "off" control |
| `rab7a_full_gtp` | nucleotide state | "on" state |
| `capza1_monomer` | monomer | |
| `capzb_monomer` | monomer | |
| `capz_heterodimer` | heterodimer | CAPZA1 + CAPZB |
| `fam129b_rab5a_gtp` | binary complex | **core hypothesis** |
| `fam129b_rab5a_gdp` | binary complex | negative control |
| `fam129b_rab7a_gtp` | binary complex | **core hypothesis** |
| `fam129b_rab7a_gdp` | binary complex | negative control |
| `fam129b_capz` | ternary complex | **core hypothesis**, h20-gpu |

Previously completed: `rab5a_full_monomer`, `rab5a_full_gdp_mg`,
`rab5a_full_gtp_mg`, `rab5a_gdp_mg_rabaptin5_complex`,
`rab5a_gtp_mg_rabaptin5_complex`, `rabep1_rabaptin5_monomer`,
`fam129b_niban2_monomer`.

---

### Stage 2 — Structural Verification
Run on the completed RAB5A-Rabaptin5 controls first. If controls pass,
apply the same analysis to FAM129B complexes.

**Why iptm alone is not enough:** iptm measures confidence, not binding
affinity. Need quantitative interface metrics to confirm AF3 is correctly
capturing GTP-state preference (real effectors bind GTP > GDP).

**Pass criterion for controls:** ΔG(GTP) < ΔG(GDP) **and** BSA(GTP) >
BSA(GDP) for RAB5A-Rabaptin5.

**Metrics computed (`src/scoring/`):**

| Metric | Tool | Script |
|--------|------|--------|
| iptm, ptm, PAE | AF3 JSON output | `parse_af3.py` |
| Interface BSA | Biopython SASA | `interface.py` |
| Residue-residue contact count | Biopython NeighborSearch (5Å cutoff) | `interface.py` |
| Switch I/II conformation | Biopython distance measurement | `interface.py` |
| PRODIGY ΔG | prodigy-prot (open-source) | `prodigy.py` |

**Output:** `results/af3_scores.csv` — all runs, all metrics in one table.

---

### Stage 3 — Interface Extraction
From the best-scoring FAM129B complex (highest iptm, confirmed by ΔG):

1. Identify contact residues on RAB5A within 5Å of FAM129B
2. Compute centroid of contact patch → X, Y, Z grid center
3. Compute extent of patch + 4Å padding → grid box dimensions
4. Strip FAM129B chain → isolated RAB5A receptor file

Both RAB5A and FAM129B will be prepared as receptors (screen from both
sides; decide which is more druggable after seeing the interface geometry).

---

### Stage 4 — Receptor Preparation
| Step | Tool | License |
|------|------|---------|
| CIF → PDB, add missing atoms/H | `pdbfixer` (OpenMM) | MIT |
| Add charges, convert to PDBQT | `meeko` (Scripps) | LGPL |

Script: `src/docking/prepare.py`

---

### Stage 5 — Ligand Library & Prep
**Library:** Enamine PPI-focused set (~5,000 compounds), free for academic
use. Downloaded automatically by the pipeline — no action needed from user.

**Custom compounds:** Add SMILES to `ligands/custom.smi` at any time; they
are included in every run automatically.

**Prep:** RDKit (BSD) — SMILES → 3D conformer (ETKDGv3) → MMFF94 minimize
→ meeko → PDBQT. Script: `src/docking/prepare.py`

---

### Stage 6 — Virtual Screening
| Tool | Role | License |
|------|------|---------|
| AutoDock Vina | Primary screen, all ~5K compounds | Apache 2.0 |
| GNINA | Rescore top 500 Vina hits | Apache 2.0 |

Script: `src/docking/screen.py`, `src/docking/gnina.py`

Both tools run against the RAB5A-side receptor. FAM129B-side receptor
run deferred until interface geometry is reviewed.

---

### Stage 7 — Triage
1. Run 200 random decoys through Vina → mean and std of score distribution
2. Convert all scores: `Z = (score − mean) / std`
3. Keep compounds with Z ≤ −2.5
4. Cross-check: compounds scoring well in both Vina **and** GNINA →
   high-confidence hits
5. Visual inspection of top 10–20 poses in notebook

Script: `src/docking/triage.py`
Output: `results/docking_hits.csv`

---

### Stage 8 — Wet-lab Handoff
Top 5–10 compounds ordered for:
- SPR or thermal shift assay — does it bind RAB5A in isolation?
- AlphaScreen or TR-FRET — does it disrupt RAB5A–FAM129B in solution?

---

## Validation & Visualization

Every stage produces plots saved to `results/plots/` and displayed inline
in `notebooks/pipeline_review.ipynb`. The notebook is the human review
checkpoint — run cells top to bottom, stop and inspect between stages.

### Stage 2 plots (`src/viz/af3_plots.py`)
| Plot | Pass signal |
|------|-------------|
| PAE heatmap per job | Low-value block at chain-chain boundary |
| pLDDT per residue | Drops flag disordered regions of FAM129B |
| GTP vs GDP ΔG bar chart | GTP bar more negative |
| GTP vs GDP BSA bar chart | GTP bar larger |
| Summary score table (all runs) | Quick cross-run sanity check |

### Stage 2–3 structure plots (`src/viz/structure_plots.py`, py3Dmol)
| Plot | What it shows |
|------|--------------|
| Full complex, chain-colored | Complex assembled correctly |
| pLDDT-colored structure | Blue=confident, red=disordered; standard AF3 QC view |
| Binding site closeup | Interface residues as sticks/spheres on transparent surface |
| GTP/GDP overlay | Switch I/II loops superimposed; shows conformational change |
| Grid box on receptor | Docking search box on the interface; confirms placement |
| Top docked pose | Best compound sitting in binding pocket |

### Stage 3 plots (`src/viz/interface_plots.py`)
| Plot | Pass signal |
|------|-------------|
| Contact residue heatmap | Dense block near Switch I/II |
| Grid box 2D diagram | Box centered on contact patch |

### Stage 6–7 plots (`src/viz/docking_plots.py`)
| Plot | What to look for |
|------|-----------------|
| Score distribution histogram | Bell curve; hits in left tail |
| Z-score ranked plot with cutoff | Clean separation of top hits |
| Vina vs GNINA scatter | Top hits cluster top-left (consensus) |
| Top 20 hits bar chart | Final shortlist with both scores |

---

## Repo Structure

```
af3-pipeline/
├── protein/
│   ├── sequences/          FASTA files (all 6 proteins)
│   ├── inputs/             AF3 JSON inputs (all 18 jobs)
│   └── templates/          JSON templates for new jobs
│
├── jobs/                   SLURM scripts (all 18 jobs)
│   └── logs/               gitignored
│
├── src/
│   ├── scoring/
│   │   ├── parse_af3.py    extract pLDDT/PAE/iptm from output JSONs
│   │   ├── interface.py    BSA, contact count, Switch I/II distances
│   │   └── prodigy.py      PRODIGY ΔG wrapper
│   ├── docking/
│   │   ├── prepare.py      CIF→PDB (pdbfixer), PDB→PDBQT (meeko)
│   │   ├── screen.py       AutoDock Vina batch screening
│   │   ├── gnina.py        GNINA rescore
│   │   └── triage.py       Z-score normalization, hit table
│   └── viz/
│       ├── af3_plots.py    PAE heatmap, pLDDT, ΔG/BSA bar charts
│       ├── structure_plots.py  3D views via py3Dmol
│       ├── interface_plots.py  contact heatmap, grid box diagram
│       └── docking_plots.py    score distributions, rankings, scatter
│
├── ligands/
│   ├── library/            PPI compound set — gitignored (large)
│   ├── custom.smi          user-provided specific compounds — tracked
│   └── decoys.smi          200 random decoys for Z-score baseline
│
├── results/
│   ├── af3_scores.csv      all AF3 runs scored — tracked
│   ├── docking_hits.csv    top compounds after triage — tracked
│   └── plots/              generated figures — tracked
│       ├── af3/
│       ├── structures/
│       ├── interface/
│       └── docking/
│
├── notebooks/
│   ├── pocket_verify.ipynb     existing
│   └── pipeline_review.ipynb  one-stop validation notebook
│
├── environment.yml         all deps: biopython, rdkit, meeko, vina,
│                           prodigy-prot, pdbfixer, py3Dmol, matplotlib
├── claude/                 planning docs (this file)
└── README.md
```

**Gitignored (large/regenerable):**
- `ligands/library/`
- `results/docking/` (raw Vina pose files)
- AF3 raw outputs live at `../output/` (outside repo) — set via env var
  `AF3_OUTPUT_DIR`, defaults to `../output` relative to repo root

---

## What the User Needs to Provide

- **Nothing required** — ligand library is downloaded automatically
- **Optional:** add specific SMILES to `ligands/custom.smi` to test
  particular compounds of interest (known drugs, published inhibitors, etc.)
- **VPN/cluster access** to submit the 11 queued SLURM jobs

---

## Build Order (when ready to proceed)

1. `environment.yml` — all dependencies pinned
2. `src/scoring/` — works on data already in hand
3. `src/viz/af3_plots.py` + `structure_plots.py` — visualize current results immediately
4. `notebooks/pipeline_review.ipynb` sections 1–3 (AF3 QC + controls)
5. `src/docking/` + remaining viz
6. `notebooks/pipeline_review.ipynb` sections 4–5 (docking + triage)
