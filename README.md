# FAM129B Screening Pipeline

**Goal:** Identify small-molecule inhibitors of the RAB5A–FAM129B
protein-protein interface using AlphaFold 3 co-folding + AutoDock Vina
virtual screening. Relevant to integrin / PD-L1 receptor trafficking in
cancer.

**Current branch:** `pipeline-build` — full pipeline code, ready to use.
**Main branch:** `main` — original manual walkthrough notes only.

---

## Project status

### FAM129B pipeline
| Stage | Status | Notes |
|-------|--------|-------|
| AF3 inputs (all proteins) | ✅ Done | All JSONs + SLURM scripts in `jobs/` |
| AF3 runs — RAB5A controls | ✅ Done | Results in `../output/` |
| AF3 runs — FAM129B complexes | ⏳ Pending | Submit jobs when on VPN |
| Environment setup | ⏳ Pending | Run `conda env create` (see below) |
| Notebook validation (Sections 0–2) | ⏳ Pending | Works with current data |
| Docking setup + screening | ⏳ Pending | Needs FAM129B complex output first |

### TRPA1 antagonist screen (new)
| Stage | Status | Notes |
|-------|--------|-------|
| Download 6V9W cryo-EM structure | ⏳ Start here — no VPN | `results/trpa1/receptor/6V9W.pdb` |
| Prep receptor PDBQT + docking grid | ⏳ After download | `src/docking/prepare.py` |
| Benchmark: HC-030031 vs 6V9W | ⏳ Validates box placement | PubChem CID 2723949 |
| Phase 1: Vina screen (~5–10K compounds) | ⏳ After grid confirmed | `src/docking/screen.py` |
| Phase 1: GNINA rescore + triage | ⏳ After Vina | `src/docking/triage.py` |
| Phase 2: AF3 co-fold top 20 hits | ⏳ Needs cluster + VPN | `jobs/trpa1_tm_hc030031.job` template |
| Phase 3: Wet-lab handoff | ⏳ After dual-validated shortlist | |

---

## 1. Environment setup

**One-time setup. Run this first before anything else.**

```bash
cd /Users/jetyue04/af3/af3-pipeline
conda env create -f environment.yml
conda activate af3-pipeline
```

This installs everything: biopython, rdkit, jupyter, py3Dmol, vina, etc.
Uses conda-forge for ARM64-native builds (Apple Silicon compatible).

To verify the install worked:

```bash
conda activate af3-pipeline
python -c "import pandas, Bio, rdkit, py3Dmol; print('OK')"
```

> **If `conda env create` fails:** run `conda env remove -n af3-pipeline`
> first to clean up, then try again.

---

## 2. Repo layout

```
af3-pipeline/
│
├── protein/
│   ├── sequences/          FASTA files for all 6 proteins
│   ├── inputs/             AF3 input JSONs (18 jobs total)
│   └── templates/          JSON templates for building new jobs
│
├── jobs/                   SLURM job scripts (one per AF3 run)
│   └── logs/               SLURM stdout/stderr (gitignored)
│
├── src/                    Pipeline code
│   ├── scoring/            Parse AF3 outputs, compute interface metrics
│   │   ├── parse_af3.py    Read summary_confidences.json → score table
│   │   ├── interface.py    BSA, contact count, Switch I/II distances
│   │   └── prodigy.py      PRODIGY binding energy (ΔG) wrapper
│   ├── docking/            Virtual screening pipeline
│   │   ├── prepare.py      CIF → PDB → PDBQT (pdbfixer + meeko)
│   │   ├── screen.py       AutoDock Vina batch screening
│   │   ├── gnina.py        GNINA deep-learning rescore
│   │   └── triage.py       Z-score normalization, hit selection
│   └── viz/                All plot generation
│       ├── af3_plots.py    PAE heatmaps, pLDDT, ΔG/BSA comparisons
│       ├── structure_plots.py  3D views via py3Dmol
│       ├── interface_plots.py  Contact maps, docking grid diagram
│       └── docking_plots.py    Score distributions, hit rankings
│
├── ligands/
│   ├── library/            Downloaded PPI compound library (gitignored)
│   ├── custom.smi          Your specific compounds to test (add here)
│   └── decoys.smi          20 decoy compounds for Z-score calibration
│
├── results/
│   ├── af3_scores.csv      Scoring results from all AF3 runs (auto-generated)
│   ├── docking_hits.csv    Top docking hits after triage (auto-generated)
│   └── plots/              All generated figures
│       ├── af3/            PAE heatmaps, pLDDT, score comparisons
│       ├── structures/     3D structure PNG exports
│       ├── interface/      Contact maps, grid box diagrams
│       └── docking/        Score distributions, hit rankings
│
├── notebooks/
│   └── pipeline_review.ipynb   Main validation notebook (start here)
│
├── claude/                 Planning docs (reference only)
│   ├── plan.md             Agreed pipeline plan — read this for full context
│   ├── guide.md            Research background and protocol
│   └── pipeline.md         Repo design decisions
│
└── environment.yml         Conda environment spec
```

---

## 3. AF3 jobs — what's been run, what's pending

### Already completed (results in `../output/`)

| Job | What it is | Key score |
|-----|-----------|-----------|
| `rab5a_full_monomer` | RAB5A alone | ptm=0.76 ✅ |
| `rab5a_full_gtp_mg` | RAB5A + GTP + Mg²⁺ | iptm=0.97 ✅ |
| `rab5a_full_gdp_mg` | RAB5A + GDP + Mg²⁺ | iptm=0.97 ✅ |
| `rab5a_gtp_mg_rabaptin5_complex` | RAB5A-GTP + Rabaptin5 (Tier 1 control) | iptm=0.70 ✅ |
| `rab5a_gdp_mg_rabaptin5_complex` | RAB5A-GDP + Rabaptin5 (negative control) | iptm=0.69 ✅ |
| `rabep1_rabaptin5_monomer` | Rabaptin5 alone | ptm=0.23 (expected — 60% disordered) |
| `fam129b_niban2_monomer` | FAM129B alone | ptm=0.75 ✅ |

**Control result:** GTP complex scores slightly better than GDP for a known
effector — AF3 is capturing nucleotide-state preference correctly.

### Pending — submit when on VPN

Connect to cluster VPN, then:

```bash
# Copy repo to cluster first (if not already there)
rsync -av /Users/jetyue04/af3/af3-pipeline/ <user>@<cluster>:/work/<user>/af3-pipeline/

# Submit all 11 jobs
cd /work/<user>/af3-pipeline
sbatch jobs/rab7a_monomer.job
sbatch jobs/rab7a_full_gdp.job
sbatch jobs/rab7a_full_gtp.job
sbatch jobs/capza1_monomer.job
sbatch jobs/capzb_monomer.job
sbatch jobs/capz_heterodimer.job
sbatch jobs/fam129b_rab5a_gtp.job      # ← core hypothesis, highest priority
sbatch jobs/fam129b_rab5a_gdp.job      # ← negative control
sbatch jobs/fam129b_rab7a_gtp.job      # ← core hypothesis
sbatch jobs/fam129b_rab7a_gdp.job      # ← negative control
sbatch jobs/fam129b_capz.job           # ← core hypothesis (uses h20-gpu)
```

Monitor jobs:
```bash
squeue -u $USER
```

When jobs finish, copy outputs back locally:
```bash
rsync -av <user>@<cluster>:/work/<user>/AF3_Output/ /Users/jetyue04/af3/output/
```

---

## 4. Running the validation notebook

The notebook is the main human-in-the-loop review step. Open it after
setting up the environment.

```bash
conda activate af3-pipeline
cd /Users/jetyue04/af3/af3-pipeline
jupyter lab notebooks/pipeline_review.ipynb
```

Or open `notebooks/pipeline_review.ipynb` directly in VS Code with the
Jupyter extension.

### What you can run right now (with existing data)

| Section | What it does | Can run now? |
|---------|-------------|--------------|
| **0 — Setup** | Imports, sets paths, lists available jobs | ✅ Yes |
| **1 — Monomer QC** | pLDDT plots + 3D views for each protein | ✅ Yes (3 of 6 monomers available) |
| **2 — Control validation** | GTP vs GDP comparison for RAB5A-Rabaptin5 | ✅ Yes |
| **3 — FAM129B complexes** | Same analysis for novel hypothesis pairs | ⏳ After SLURM jobs |
| **4 — Docking setup** | Interface extraction, receptor prep, grid box | ⏳ After Section 3 |
| **5 — Screening results** | Score distributions, top hits, 3D poses | ⏳ After docking run |

### How to step through it

Run cells top to bottom. Each section ends with a **PASS/FAIL verdict**
printed in the cell output:
- `✓ PASS` — safe to continue
- `✗ FAIL` or `⚠ WARNING` — read the message before proceeding

**Do not interpret FAM129B results if Section 2 fails.** The controls must
pass first.

---

## 5. Adding your own compounds to screen

Edit `ligands/custom.smi` — one compound per line, format `SMILES name`:

```
CC1=CC=C(C=C1)NC(=O)C compound_001
c1ccc2c(c1)cccc2 compound_002
```

These are automatically included in every docking run alongside the
PPI library. You can add compounds at any time, even after the pipeline
is already set up.

---

## 6. Full pipeline walkthrough (end to end)

Once all AF3 jobs are back and the notebook controls pass:

### Step 1 — Score all AF3 outputs

```bash
conda activate af3-pipeline
cd /Users/jetyue04/af3/af3-pipeline
python src/scoring/parse_af3.py
# → writes results/af3_scores.csv
```

### Step 2 — Open notebook, run Sections 0–3

Inspect plots. Identify the best-scoring FAM129B complex (highest iptm,
GTP > GDP confirmed). Note the job name — you'll need it for docking.

### Step 3 — Prepare receptor for docking

```bash
# Replace JOB_NAME with the best FAM129B complex job
python -c "
from src.docking.prepare import prepare_receptor
prepare_receptor(
    '../output/JOB_NAME/JOB_NAME_model.cif',
    chain_id='B',        # RAB5A chain
    output_dir='results/docking'
)
"
# → writes results/docking/chain_B.pdbqt
```

### Step 4 — Run notebook Section 4

Verify the docking grid box is centered on the binding interface before
screening. Do not skip this — a misplaced grid box wastes the whole run.

### Step 5 — Screen the library

```bash
python src/docking/screen.py
# Screens ligands/library/ against the receptor
# → writes results/docking_hits.csv
# Takes ~1–4 hours depending on library size
```

### Step 6 — Run notebook Section 5

Review score distribution, Z-score rankings, top hits. Inspect 3D poses
of the top compounds. Shortlist 5–10 for wet-lab ordering.

---

## 7. Scores — what they mean

| Metric | What it measures | Good value |
|--------|-----------------|------------|
| `ptm` | Overall structure confidence (monomer) | > 0.7 |
| `iptm` | Interface confidence (complex) | > 0.75 = strong, 0.5–0.75 = weak |
| `fraction_disordered` | % of residues with low pLDDT | < 0.3 ideally |
| `bsa` | Buried surface area at interface (Å²) | > 500 Å² suggests real contact |
| `prodigy_dg` | Predicted binding energy (kcal/mol) | More negative = tighter binding |
| `vina_zscore` | Z-score vs random decoys | ≤ −2.5 = hit |
| `consensus_flag` | Hit in both Vina AND GNINA | True = high confidence |

---

## 8. Reference

Full pipeline plan with scientific rationale: [`claude/plan.md`](claude/plan.md)

Research background (proteins, biology, protocol): [`claude/guide.md`](claude/guide.md)
