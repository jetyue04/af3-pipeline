# TRPA1 Antagonist Screen — Plan

**Goal:** Identify novel small-molecule antagonists of human TRPA1 (O75762) using
the experimental cryo-EM structure + AutoDock Vina virtual screening.
AF3 co-folding is used at the end to validate top hits, not for the primary screen.

**Why this target:** TRPA1 is a non-selective cation channel gated by noxious stimuli
(cold, reactive chemicals, oxidative stress). It drives pain and neurogenic
inflammation. A selective small-molecule antagonist has therapeutic value for pain,
asthma, and itch.

**Why use the PDB structure instead of AF3 for screening:**
- TRPA1 has high-resolution experimental cryo-EM structures with ligands already bound
- An experimental structure with a co-crystallized ligand (6V9W) directly defines the
  pocket geometry — more reliable for docking than a de-novo prediction
- No cluster time or VPN required for the primary screen
- AF3 co-folding is reserved for top-hit validation, where per-compound compute is justified

---

## Phase 0 — Receptor preparation (start here)

### Structure to use: PDB 6V9W
- Human TRPA1, 3.0 Å resolution cryo-EM
- Co-crystallized with **A-967079** (a known selective TRPA1 antagonist)
- The ligand defines the exact pocket geometry — ideal for docking grid placement
- Chain A is the subunit of interest

### Alternative: PDB 7UM5
- Human TRPA1 apo, 2.7 Å — use if you want a ligand-unbiased receptor conformation
- Higher resolution but no ligand to anchor the docking box

### Steps

**Step 1 — Download the structure**
```bash
cd /Users/jetyue04/af3/af3-pipeline
mkdir -p results/trpa1/receptor

# Download 6V9W from RCSB
curl -o results/trpa1/receptor/6V9W.pdb \
  "https://files.rcsb.org/download/6V9W.pdb"
```

**Step 2 — Prepare receptor PDBQT**

Use the existing `src/docking/prepare.py` — same script as FAM129B:
```bash
conda activate af3-pipeline

python -c "
from src.docking.prepare import prepare_receptor
prepare_receptor(
    'results/trpa1/receptor/6V9W.pdb',
    chain_id='A',                        # one TRPA1 subunit
    output_dir='results/trpa1/receptor'
)
"
# → writes results/trpa1/receptor/chain_A.pdbqt
```

**Step 3 — Define the docking grid box**

The A-967079 ligand in 6V9W sits in the TM4–S5 linker / S6 intracellular pocket.
Extract the ligand coordinates to center the grid:

```bash
python -c "
from Bio.PDB import PDBParser
import numpy as np

parser = PDBParser(QUIET=True)
structure = parser.get_structure('6V9W', 'results/trpa1/receptor/6V9W.pdb')
# A-967079 is residue name '967' or similar — check with grep first
import subprocess
subprocess.run(['grep', 'HETATM', 'results/trpa1/receptor/6V9W.pdb'])
"
```

Note the ligand residue name from the output, then compute centroid:
```bash
python -c "
from Bio.PDB import PDBParser
import numpy as np

parser = PDBParser(QUIET=True)
s = parser.get_structure('6V9W', 'results/trpa1/receptor/6V9W.pdb')
# Replace 'LIG' below with the actual residue name of A-967079 in 6V9W
lig_atoms = [a for a in s.get_atoms() if a.parent.resname == 'LIG']
coords = np.array([a.coord for a in lig_atoms])
center = coords.mean(axis=0)
extent = coords.max(axis=0) - coords.min(axis=0)
box = extent + 8   # 4 Å padding each side
print(f'Grid center: x={center[0]:.1f}  y={center[1]:.1f}  z={center[2]:.1f}')
print(f'Grid size:   x={box[0]:.1f}   y={box[1]:.1f}   z={box[2]:.1f}')
"
```

Record these values — they go into `src/docking/screen.py` as `center_x/y/z` and `size_x/y/z`.

---

## Phase 1 — Ligand library screen (AutoDock Vina)

### Library
**Enamine ion-channel focused set** (~5–10K compounds) — better chemistry for a TRP
channel than the PPI library used for FAM129B.

Download via Enamine's academic access portal, or reuse the existing PPI library in
`ligands/library/` as a fast first pass (less optimal but immediately available).

Also add any specific SMILES of interest to `ligands/custom.smi`.

### Screen
```bash
conda activate af3-pipeline
python src/docking/screen.py \
    --receptor results/trpa1/receptor/chain_A.pdbqt \
    --ligand_dir ligands/library/ \
    --output_dir results/trpa1/docking \
    --center_x X --center_y Y --center_z Z \
    --size_x SX --size_y SY --size_z SZ
# → writes results/trpa1/docking_hits.csv
# Runtime: 1–4 hrs depending on library size
```

*(Replace X/Y/Z and SX/SY/SZ with the values from Phase 0 Step 3.)*

### GNINA rescore (top 500 Vina hits)
```bash
python src/docking/gnina.py \
    --receptor results/trpa1/receptor/chain_A.pdbqt \
    --poses_dir results/trpa1/docking \
    --output results/trpa1/gnina_scores.csv
```

### Triage
```bash
python src/docking/triage.py \
    --scores results/trpa1/docking_hits.csv \
    --gnina results/trpa1/gnina_scores.csv \
    --output results/trpa1/top_hits.csv
```

**Hit criteria** (same thresholds as FAM129B pipeline):
| Metric | Threshold |
|--------|-----------|
| `vina_zscore` | ≤ −2.5 |
| `consensus_flag` | True (hit in both Vina AND GNINA) |

---

## Phase 2 — AF3 co-fold validation (top 20 hits only)

After Phase 1 shortlists ~20 candidates, run each as an AF3 co-fold job with the
TRPA1 TM domain to get a high-confidence second opinion on the binding mode.

### Construct for co-folding
- **TRPA1 TM domain: residues 682–1119** (438 aa)
  - Starts after the last ankyrin repeat, just before S1 helix
  - Covers S1–S6, TRP helix, C-terminal coiled-coil
  - Captures the full antagonist binding pocket at 40% of full-length
- Fetch and trim sequence:

```bash
curl "https://rest.uniprot.org/uniprotkb/O75762.fasta" \
  > protein/sequences/O75762_TRPA1.fasta

python3 -c "
from Bio import SeqIO
record = next(SeqIO.parse('protein/sequences/O75762_TRPA1.fasta', 'fasta'))
tm_seq = str(record.seq[681:1119])   # 0-indexed: residues 682-1119
print(f'TM domain: {len(tm_seq)} aa (expected 438)')
print(tm_seq)
"
```

### JSON format for a hit compound (one file per compound)
Each hit gets its own JSON — copy `protein/inputs/trpa1_tm_hc030031.json` as a
template and replace the SMILES:

```json
{
  "name": "TRPA1_TM_682_1119_hit_001",
  "modelSeeds": [1, 2, 3],
  "sequences": [
    {
      "protein": {
        "id": "A",
        "sequence": "PASTE_TM_SEQUENCE_HERE"
      }
    },
    {
      "ligand": {
        "id": "B",
        "smiles": "SMILES_OF_HIT_COMPOUND"
      }
    }
  ],
  "dialect": "alphafold3",
  "version": 2
}
```

Three model seeds (`[1, 2, 3]`) per compound — consistent placement across
seeds = confident pose.

### Co-fold pass criteria
| Metric | Threshold | Meaning |
|--------|-----------|---------|
| `iptm` | > 0.5 | Ligand–protein interface confidence |
| Ligand PAE vs TM helices | < 10 Å | Pose is spatially confident |
| Visual position | Between TM helices, intracellular face | Not floating outside |
| Consistent across 3 seeds | ≥ 2/3 agree | Reproducible binding mode |

Hits that pass both Phase 1 (Vina/GNINA) AND Phase 2 (AF3 iptm > 0.5 in ≥ 2/3 seeds)
are the **high-confidence shortlist** for wet-lab ordering.

---

## Phase 3 — Wet-lab validation handoff

Top 5–10 dual-validated compounds ordered for:
- **Thermal shift assay** — does it bind TRPA1 in isolation?
- **FLIPR calcium flux assay** — does it block TRPA1 channel activity in cells?
- **Whole-cell patch clamp** (optional) — direct electrophysiology confirmation

---

## Benchmark run: HC-030031 (do this before Phase 1)

Before screening unknowns, run one known antagonist through the full pipeline as a
positive control. If HC-030031 scores as a top hit against 6V9W, the docking box
and pipeline are correctly configured.

- HC-030031 SMILES: `CC(C)c1ccc(NC(=O)Cn2cnc3c(=O)n(C)c(=O)n(C)c32)cc1`
  - **Verify from PubChem CID 2723949 before use**
- Add it to `ligands/custom.smi` with a name so it appears in the triage table
- Expected result: Vina Z-score ≤ −2.5, high GNINA score

---

## Files created

| File | Purpose |
|------|---------|
| `protein/inputs/trpa1_monomer.json` | AF3 monomer QC (optional, parallel) |
| `protein/inputs/trpa1_tm_hc030031.json` | AF3 co-fold template for hit validation |
| `jobs/trpa1_monomer.job` | SLURM script for monomer QC |
| `jobs/trpa1_tm_hc030031.job` | SLURM script for co-fold validation |

---

## Status

| Step | Status |
|------|--------|
| Download 6V9W from RCSB | ⏳ Do this first — no VPN needed |
| Prep receptor PDBQT (chain A) | ⏳ After download |
| Define docking grid box | ⏳ After receptor prep |
| Benchmark run: HC-030031 | ⏳ Validates box placement |
| Phase 1: Vina screen (~5–10K compounds) | ⏳ After box confirmed |
| Phase 1: GNINA rescore + triage | ⏳ After Vina |
| Phase 2: AF3 co-fold top 20 hits | ⏳ After triage; needs cluster + VPN |
| Phase 3: Wet-lab handoff | ⏳ After dual-validated shortlist |
