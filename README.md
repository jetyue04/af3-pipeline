# AF3 Virtual Screening Pipeline

Two active targets sharing the same infrastructure:

| Target | Type | Goal |
|--------|------|------|
| **FAM129B × RAB5A** | Protein–protein interaction | Inhibitors of FAM129B–RAB5A interface; integrin/PD-L1 trafficking in cancer |
| **TRPA1** | Ion channel | Novel antagonists binding the TM-domain antagonist pocket |

**Active branch:** `trpa1-docking` — TRPA1 screen in progress.
**Stable branch:** `main` — FAM129B AF3 analysis complete.

---

## Project status

### FAM129B × RAB5A

| Stage | Status | Notes |
|-------|--------|-------|
| AF3 controls (RAB5A monomers, Rabaptin5) | ✅ Done | Results in `../output/` |
| AF3 — full FAM129B + RAB5A GTP/GDP | ✅ Done | No interface (disordered N-term masked signal) |
| AF3 — HB domain (192–576) + RAB5A GTP | ✅ Done | **pair_iptm=0.56** — weak/candidate interface |
| AF3 — HB domain + RAB5A GDP | ✅ Done | pair_iptm=0.23 — confirms GTP-selectivity |
| AF3 — Motif subdomains (Motif1, Motif2) × RAB5A | ✅ Done | Both fail; cooperative binding required |
| Summary report | ✅ Done | `results/fam129b_rab5a_summary_report.md` |
| RAB7A, CapZ complex jobs | ⏳ Submit when on VPN | Jobs ready in `jobs/` |
| Docking screen (PPI inhibitors) | ⏳ After RAB5A/RAB7A complex outputs | |

### TRPA1 antagonist screen

| Stage | Status | Notes |
|-------|--------|-------|
| Download PDB 6V9W (cryo-EM + A-967079) | ⏳ Start here — no VPN | `results/trpa1/receptor/6V9W.pdb` |
| Prepare chain A as PDBQT | ⏳ After download | `src/docking/prepare.py` |
| Define docking grid box | ⏳ After receptor prep | Derived from A-967079 ligand coords |
| Benchmark: HC-030031 (validates box) | ⏳ After grid defined | Expected vina_zscore ≤ −2.5 |
| Phase 1: Vina screen (~5–10K compounds) | ⏳ After benchmark passes | `src/docking/screen.py` |
| Phase 1: GNINA rescore + triage | ⏳ After Vina | `src/docking/gnina.py`, `triage.py` |
| Phase 2: AF3 co-fold top 20 hits | ⏳ Needs cluster + VPN | TM domain 682–1119 + ligand SMILES |
| Phase 3: Wet-lab handoff | ⏳ After dual-validated shortlist | Thermal shift + FLIPR assay |

---

## 1. Environment setup

**One-time setup. Run this first.**

```bash
cd /Users/jetyue04/af3/af3-pipeline
conda env create -f environment.yml
conda activate af3-pipeline
```

Installs: biopython, rdkit, jupyter, py3Dmol, vina, gnina, meeko, pdbfixer.
Uses conda-forge ARM64-native builds (Apple Silicon compatible).

Verify:
```bash
conda activate af3-pipeline
python -c "import pandas, Bio, rdkit, py3Dmol; print('OK')"
```

> If `conda env create` fails: `conda env remove -n af3-pipeline` then retry.

---

## 2. Repo layout

```
af3-pipeline/
│
├── protein/
│   ├── sequences/          FASTA files (FAM129B, RAB5A, RAB7A, CAPZA1, CAPZB, TRPA1)
│   └── inputs/             AF3 input JSONs (22 jobs total)
│
├── jobs/                   SLURM job scripts (one per AF3 run)
│   └── logs/               SLURM stdout/stderr (gitignored)
│
├── src/
│   ├── scoring/
│   │   ├── parse_af3.py    Read summary_confidences.json → score table
│   │   ├── interface.py    BSA, contact count, switch region distances
│   │   └── prodigy.py      PRODIGY binding energy wrapper
│   ├── docking/
│   │   ├── prepare.py      CIF/PDB → PDBQT (pdbfixer + meeko)
│   │   ├── screen.py       AutoDock Vina batch screen
│   │   ├── gnina.py        GNINA deep-learning rescore
│   │   └── triage.py       Z-score normalization, hit selection
│   └── viz/
│       ├── af3_plots.py        PAE heatmaps, pLDDT, score comparisons
│       ├── structure_plots.py  3D views via py3Dmol
│       ├── interface_plots.py  Contact maps, docking grid diagrams
│       └── docking_plots.py    Score distributions, hit rankings
│
├── ligands/
│   ├── library/            Compound library for screening (gitignored — download separately)
│   ├── custom.smi          Add your own SMILES here (one per line: SMILES name)
│   └── decoys.smi          20 decoys for Z-score calibration
│
├── results/
│   ├── af3_scores.csv                      All AF3 run scores
│   ├── fam129b_rab5a_summary_report.md     PI-ready FAM129B results report
│   ├── view_hb_gtp.cxc                     ChimeraX: HB domain + RAB5A GTP
│   ├── view_motifs_comparison.cxc          ChimeraX: motif subdomains vs HB domain
│   └── trpa1/                              TRPA1 screen outputs (created by notebook)
│       ├── receptor/       6V9W.pdb, chain_A.pdbqt
│       ├── docking/        Per-compound pose files + scores
│       └── af3_validation/ Top-hit AF3 co-fold results
│
├── notebooks/
│   ├── pipeline_review.ipynb   FAM129B AF3 result review
│   └── trpa1_screen.ipynb      TRPA1 docking screen — step-by-step workflow
│
├── claude/                 Planning docs (reference)
│   ├── plan.md             FAM129B pipeline plan
│   ├── trpa1_plan.md       TRPA1 screen plan (detailed)
│   └── guide.md            Research background
│
└── environment.yml
```

---

## 3. FAM129B × RAB5A — key results

Full report: [`results/fam129b_rab5a_summary_report.md`](results/fam129b_rab5a_summary_report.md)

### AF3 scores (all completed runs)

| Construct | pair_iptm | PAE_min | Verdict |
|-----------|-----------|---------|---------|
| Full FAM129B + RAB5A GDP | 0.28 | 19.9 Å | No interface |
| Full FAM129B + RAB5A GTP | 0.16 | 24.8 Å | No interface |
| HB domain (192–576) + RAB5A GDP | 0.23 | 15.1 Å | No interface |
| **HB domain + RAB5A GTP** | **0.56** | **6.2 Å** | **Candidate interface ✓** |
| Motif 1 (216–287) + RAB5A GTP | 0.23 | — | Disordered alone |
| Motif 2 (387–448) + RAB5A GTP | 0.19 | — | Disordered alone |

### Key conclusions

1. **HB domain (192–576) is the minimal binding unit.** Removing the disordered N-terminus (1–191) was critical — it masked the binding signal in the full-length runs.
2. **GTP-selective interaction.** pair_iptm jumps from 0.23 (GDP) to 0.56 (GTP) for the HB domain — a 2.4× increase driven by Switch II accessibility.
3. **Cooperative/bivalent binding.** Two contact patches (226–277 and 397–438) are both required simultaneously. Neither binds alone as an isolated peptide.
4. **Critical residue: FAM129B Q403.** Only residue contacting RAB5A Switch II (E72, 4.4 Å). The GTP-selectivity sensor.
5. **Recommended first experiment:** Q403A mutagenesis + co-IP with RAB5A Q79L (constitutively GTP-loaded).

### Pending jobs (submit on VPN)

```bash
# rsync repo to cluster first, then:
sbatch jobs/fam129b_rab7a_gtp.job
sbatch jobs/fam129b_rab7a_gdp.job
sbatch jobs/fam129b_capz.job
sbatch jobs/rab7a_full_gtp.job
sbatch jobs/rab7a_full_gdp.job
sbatch jobs/rab7a_monomer.job
sbatch jobs/capza1_monomer.job
sbatch jobs/capzb_monomer.job
sbatch jobs/capz_heterodimer.job
```

---

## 4. TRPA1 antagonist screen — full workflow

**Notebook:** `notebooks/trpa1_screen.ipynb` — run this cell by cell.

The screen uses an **experimental cryo-EM structure** (PDB 6V9W) rather than AF3 for the primary docking, because a co-crystallized antagonist (A-967079) already defines the exact binding pocket geometry. AF3 co-folding is reserved for hit validation only.

### Why this receptor and pocket

- **PDB 6V9W:** Human TRPA1 at 3.0 Å, co-crystallized with **A-967079** (selective antagonist)
- The A-967079 ligand sits in the **TM4–S5 linker / intracellular S6 pocket** — the canonical small-molecule antagonist site
- The ligand coordinates directly define the docking box center and size
- Alternative: PDB 7UM5 (apo, 2.7 Å) — higher resolution but no ligand anchor; use only if you suspect biased pocket geometry

### Funnel logic

```
~5–10K compounds
    → Vina screen (fast, ~seconds/compound)      [vina_zscore ≤ −2.5]
        → top 500 → GNINA rescore (accurate CNN)  [consensus_flag = True]
            → top ~20 → AF3 co-fold validation    [iptm > 0.5, ≥ 2/3 seeds]
                → 5–10 compound shortlist for wet lab
```

Each stage filters by orthogonal criteria; the final shortlist has three independent methods in agreement.

### Phase 0 — Receptor preparation

**Step 1 — Download 6V9W**
```bash
mkdir -p results/trpa1/receptor
curl -o results/trpa1/receptor/6V9W.pdb "https://files.rcsb.org/download/6V9W.pdb"
```

**Step 2 — Prepare PDBQT (chain A)**

TRPA1 is a homotetramer; dock to one subunit:
```bash
conda activate af3-pipeline
python -c "
from src.docking.prepare import prepare_receptor
prepare_receptor('results/trpa1/receptor/6V9W.pdb',
                 chain_id='A',
                 output_dir='results/trpa1/receptor')
"
# → results/trpa1/receptor/chain_A.pdbqt
```

**Step 3 — Extract docking grid box from A-967079**

Find the ligand residue name:
```bash
grep HETATM results/trpa1/receptor/6V9W.pdb | awk '{print $4}' | sort -u
```

Then compute centroid and box size (notebook Section 3 does this interactively with a 3D preview).

**Step 4 — Benchmark HC-030031**

Before screening unknowns, dock HC-030031 (known antagonist, PubChem CID 2723949) to confirm the box is correct. Expected: vina_zscore ≤ −2.5. If it fails, the box needs adjustment before proceeding.

HC-030031 SMILES (verify from PubChem before use):
```
CC(C)c1ccc(NC(=O)Cn2cnc3c(=O)n(C)c(=O)n(C)c32)cc1
```

### Phase 1 — Vina screen + GNINA rescore

```bash
conda activate af3-pipeline

# Vina screen (replace X Y Z SX SY SZ with values from Step 3)
python src/docking/screen.py \
    --receptor results/trpa1/receptor/chain_A.pdbqt \
    --ligand_dir ligands/library/ \
    --output_dir results/trpa1/docking \
    --center_x X --center_y Y --center_z Z \
    --size_x SX --size_y SY --size_z SZ
# → results/trpa1/docking_hits.csv  (~1–4 hrs)

# GNINA rescore top 500
python src/docking/gnina.py \
    --receptor results/trpa1/receptor/chain_A.pdbqt \
    --poses_dir results/trpa1/docking \
    --output results/trpa1/gnina_scores.csv

# Triage: consensus hits only
python src/docking/triage.py \
    --scores results/trpa1/docking_hits.csv \
    --gnina results/trpa1/gnina_scores.csv \
    --output results/trpa1/top_hits.csv
```

Hit criteria:

| Metric | Threshold |
|--------|-----------|
| `vina_zscore` | ≤ −2.5 |
| `consensus_flag` | True (Vina AND GNINA) |

### Phase 2 — AF3 co-fold validation (cluster, needs VPN)

Top 20 hits from Phase 1 each get an AF3 job: TRPA1 TM domain (682–1119, 438 aa) co-folded with the hit compound as a SMILES ligand. Template: `protein/inputs/trpa1_tm_hc030031.json` — swap the SMILES, keep everything else.

The notebook (Section 8) auto-generates all 20 JSONs and SLURM job scripts from `top_hits.csv`.

Co-fold pass criteria:

| Metric | Threshold | Meaning |
|--------|-----------|---------|
| `iptm` | > 0.5 | Ligand–protein interface confidence |
| Ligand PAE vs TM helices | < 10 Å | Spatially confident pose |
| Seeds in agreement | ≥ 2 of 3 | Reproducible binding mode |

### Phase 3 — Wet-lab handoff

Top 5–10 dual-validated compounds ordered for:
- **Thermal shift assay** — direct binding confirmation
- **FLIPR calcium flux assay** — functional channel inhibition in cells
- **Whole-cell patch clamp** (optional) — gold standard electrophysiology

---

## 5. Notebooks

| Notebook | Purpose | Can run now? |
|----------|---------|--------------|
| `pipeline_review.ipynb` | FAM129B AF3 result scoring, PAE heatmaps, pLDDT | ✅ Yes (with existing output data) |
| `trpa1_screen.ipynb` | Full TRPA1 docking screen — Phase 0 → 2 | ✅ Start Phase 0 now (no VPN) |

```bash
conda activate af3-pipeline
jupyter lab
```

---

## 6. Adding your own compounds

Edit `ligands/custom.smi` — one compound per line, `SMILES name`:

```
CC1=CC=C(C=C1)NC(=O)C compound_001
CC(C)c1ccc(NC(=O)Cn2cnc3c(=O)n(C)c(=O)n(C)c32)cc1 HC-030031
```

These are included automatically in every docking run alongside the main library.

---

## 7. Scores — what they mean

| Metric | What it measures | Good threshold |
|--------|-----------------|----------------|
| `ptm` | Overall fold confidence (monomer) | > 0.7 |
| `iptm` | Interface confidence (complex) | > 0.75 strong, 0.5–0.75 candidate |
| `pair_iptm` | Per-chain-pair interface confidence | > 0.5 candidate |
| `fraction_disordered` | Residues with low pLDDT | < 0.3 |
| `bsa` | Buried surface area at interface (Å²) | > 500 Å² |
| `prodigy_dg` | Predicted binding energy (kcal/mol) | More negative = tighter |
| `vina_zscore` | Docking score vs library background | ≤ −2.5 = hit |
| `consensus_flag` | Hit in both Vina AND GNINA | True = high confidence |

---

## 8. Reference

- FAM129B pipeline plan: [`claude/plan.md`](claude/plan.md)
- TRPA1 screen plan (detailed): [`claude/trpa1_plan.md`](claude/trpa1_plan.md)
- Research background: [`claude/guide.md`](claude/guide.md)
- FAM129B results report: [`results/fam129b_rab5a_summary_report.md`](results/fam129b_rab5a_summary_report.md)
