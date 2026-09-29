# PI Guide: FAM129B Screening Pipeline
### How to run everything, step by step — no prior coding experience required

---

## Table of Contents

1. [What this project does](#1-what-this-project-does)
2. [What you will need](#2-what-you-will-need)
3. [Git basics — getting and saving the code](#3-git-basics--getting-and-saving-the-code)
4. [Connecting to the Duke cluster](#4-connecting-to-the-duke-cluster)
5. [Copying the project to the cluster](#5-copying-the-project-to-the-cluster)
6. [Understanding the JSON input files](#6-understanding-the-json-input-files)
7. [Submitting AlphaFold 3 jobs to Slurm](#7-submitting-alphafold-3-jobs-to-slurm)
8. [Monitoring your jobs](#8-monitoring-your-jobs)
9. [Copying results back to your laptop](#9-copying-results-back-to-your-laptop)
10. [Running the analysis pipeline](#10-running-the-analysis-pipeline)
11. [What the scores mean](#11-what-the-scores-mean)
12. [Full pipeline at a glance](#12-full-pipeline-at-a-glance)
13. [Troubleshooting](#13-troubleshooting)

---

## 1. What this project does

This pipeline answers one question: **which small molecules might block the interaction between the proteins FAM129B and RAB5A?**

It does this in three stages:

| Stage | Tool | What happens |
|-------|------|--------------|
| **Structural prediction** | AlphaFold 3 (AF3) | Predicts 3D structures of the proteins docked together; runs on the Duke GPU cluster |
| **Scoring & validation** | Python scripts | Measures confidence in the predicted structures; checks that GTP-state complexes score better than GDP-state (a known biology control) |
| **Virtual drug screening** | AutoDock Vina + GNINA | Tests ~5,000 small molecules against the interface; ranks and flags the best hits for wet-lab follow-up |

Everything is already set up. You mainly need to (a) submit the AF3 jobs to the cluster and (b) run a Jupyter notebook for analysis.

---

## 2. What you will need

- A **Mac** (all commands below are for macOS Terminal)
- **Duke VPN** — you must be connected whenever you access the cluster
- A **Duke NetID** — your login username for the cluster
- The project folder already on your laptop: `/Users/jetyue04/af3/af3-pipeline/`

> **How to open Terminal:** Press `Command + Space`, type `Terminal`, press Enter.

---

## 3. Git basics — getting and saving the code

Git is a version-control system. Think of it like Google Docs version history for code — every change is tracked, and you can always go back. You do **not** need to understand Git deeply to run this pipeline, but here are the three commands you will actually use.

### 3.1 Checking the current state

```bash
cd /Users/jetyue04/af3/af3-pipeline
git status
```

This shows which files have been changed. A clean output looks like:

```
On branch pipeline-build
nothing to commit, working tree clean
```

### 3.2 Downloading the latest updates from GitHub

If your student sends you an update:

```bash
cd /Users/jetyue04/af3/af3-pipeline
git pull
```

That's it. Git merges any new code automatically.

### 3.3 Saving your own changes (committing)

If you edit a file (for example, you add a compound to `ligands/custom.smi`), you can save that change permanently:

```bash
cd /Users/jetyue04/af3/af3-pipeline
git add ligands/custom.smi
git commit -m "Add compound X to screening list"
```

The message in quotes is just a short note to yourself — write anything descriptive.

### 3.4 Uploading your changes to GitHub

```bash
git push
```

### Key concept: branches

The project has two branches:
- **`pipeline-build`** — the active working branch with all the code. Use this.
- **`main`** — original notes only, ignore it.

To confirm you are on the right branch:

```bash
git branch
```

You should see `* pipeline-build` with an asterisk. If you are on a different branch:

```bash
git checkout pipeline-build
```

---

## 4. Connecting to the Duke cluster

The Duke Research Computing cluster is where AlphaFold 3 actually runs (it requires powerful GPUs that are not on your laptop). You connect to it via SSH — a secure terminal connection.

### 4.1 Connect to Duke VPN first

You **must** be on Duke VPN before connecting to the cluster. Use the Cisco AnyConnect app on your Mac and connect to `vpn.duke.edu`.

### 4.2 Open a terminal and connect

Replace `NETID` with your Duke NetID (e.g., `abc12`):

```bash
ssh NETID@dkucc-login-01.rc.duke.edu
```

Type your Duke password when prompted. You will land at a prompt that looks like:

```
[NETID@dkucc-login-01 ~]$
```

You are now on the cluster. **Do not run computational jobs directly here** — the login node is only for submitting jobs and moving files. All heavy computation goes through Slurm (the job scheduler).

### 4.3 To disconnect when done

```bash
exit
```

---

## 5. Copying the project to the cluster

The first time you use the cluster, you need to copy the pipeline code from your laptop. Run this command **on your laptop** (not on the cluster):

Replace `NETID` with your actual Duke NetID.

```bash
rsync -av /Users/jetyue04/af3/af3-pipeline/ NETID@dkucc-login-01.rc.duke.edu:/work/NETID/af3-pipeline/
```

`rsync` copies only new or changed files, so re-running this later (after a `git pull`) will sync any updates quickly.

> **What the flags mean:** `-a` preserves file permissions and timestamps; `-v` prints each file name as it copies so you can see progress.

---

## 6. Understanding the JSON input files

AlphaFold 3 needs to know **which proteins (and ligands) to fold together**. You tell it this via a JSON file — a structured text file. All 18 JSON files are pre-written and ready to use in `protein/inputs/`.

### 6.1 What a JSON file looks like

Here is the file for the core hypothesis job (`fam129b_rab5a_gtp.json`):

```json
{
  "name": "FAM129B_RAB5A_GTP_Mg_complex",
  "modelSeeds": [1],
  "sequences": [
    {
      "protein": {
        "id": "A",
        "sequence": "MGDVLSTHLD..."
      }
    },
    {
      "protein": {
        "id": "B",
        "sequence": "MASRGATRPN..."
      }
    },
    {
      "ligand": {
        "id": "C",
        "ccdCodes": ["GTP"]
      }
    },
    {
      "ligand": {
        "id": "D",
        "ccdCodes": ["MG"]
      }
    }
  ],
  "dialect": "alphafold3",
  "version": 2
}
```

### 6.2 What each field means

| Field | Meaning |
|-------|---------|
| `"name"` | A label for this run — used to name the output folder |
| `"modelSeeds"` | `[1]` means one structural prediction. Add more (e.g. `[1, 2, 3]`) to get multiple independent predictions |
| `"sequences"` | List of molecules: proteins get full amino acid sequences; ligands get their PDB chemical code |
| `"id"` | Chain letter — A, B, C, D, etc. Just a label |
| `"sequence"` | Full amino acid sequence (one-letter codes) from UniProt |
| `"ccdCodes"` | Standard chemical code: `"GTP"`, `"GDP"`, `"MG"` (magnesium), etc. |
| `"dialect"` / `"version"` | Always `"alphafold3"` and `2` — the cluster requires version 2 |

### 6.3 Pre-written JSON files (you do not need to modify these)

All 18 jobs are already prepared:

| File | What it folds | Purpose |
|------|--------------|---------|
| `rab5a_monomer.json` | RAB5A alone | Quality check |
| `rab5a_full_gtp.json` | RAB5A + GTP + Mg²⁺ | "Active" state baseline |
| `rab5a_full_gdp.json` | RAB5A + GDP + Mg²⁺ | "Inactive" state baseline |
| `fam129b_monomer.json` | FAM129B alone | Quality check |
| `fam129b_rab5a_gtp.json` | **FAM129B + RAB5A + GTP + Mg²⁺** | **Core hypothesis — active state** |
| `fam129b_rab5a_gdp.json` | FAM129B + RAB5A + GDP + Mg²⁺ | Negative control — should score lower |
| `fam129b_rab7a_gtp.json` | FAM129B + RAB7A + GTP + Mg²⁺ | Second hypothesis |
| `fam129b_rab7a_gdp.json` | FAM129B + RAB7A + GDP + Mg²⁺ | Negative control |
| `fam129b_capz.json` | FAM129B + CapZ heterodimer | Third hypothesis |
| ...and 9 more | Various controls and monomers | See `jobs/README.md` |

### 6.4 If you need to create a new JSON for a new protein

Copy the closest existing template from `protein/templates/`, then change:
1. `"name"` — give it a unique label (no spaces; use underscores)
2. `"sequence"` — paste the amino acid sequence from UniProt
3. Add or remove `"ligand"` blocks as needed

> **Important:** Do NOT add a `"description"` field or change `"version"` to 4. The Duke cluster's AF3 install does not support these and the job will fail silently.

---

## 7. Submitting AlphaFold 3 jobs to Slurm

Slurm is the job scheduler on the cluster. You tell it "run this script," and it puts it in a queue behind other users' jobs and starts it when a GPU is free.

### 7.1 SSH to the cluster and navigate to the project

```bash
ssh NETID@dkucc-login-01.rc.duke.edu
cd /work/NETID/af3-pipeline
```

### 7.2 Submit a single job

```bash
sbatch jobs/fam129b_rab5a_gtp.job
```

Slurm responds immediately with a job ID number:

```
Submitted batch job 1234567
```

Write this number down — you can use it to check status later.

### 7.3 Submit all pending jobs at once

To submit all 11 pending jobs in one go:

```bash
cd /work/NETID/af3-pipeline

sbatch jobs/rab7a_monomer.job
sbatch jobs/rab7a_full_gdp.job
sbatch jobs/rab7a_full_gtp.job
sbatch jobs/capza1_monomer.job
sbatch jobs/capzb_monomer.job
sbatch jobs/capz_heterodimer.job
sbatch jobs/fam129b_rab5a_gtp.job
sbatch jobs/fam129b_rab5a_gdp.job
sbatch jobs/fam129b_rab7a_gtp.job
sbatch jobs/fam129b_rab7a_gdp.job
sbatch jobs/fam129b_capz.job
```

> **GPU note:** Most jobs use `l20-gpu` (48 GB). The `fam129b_capz.job` is set to `h20-gpu` (96 GB) because it is a three-protein complex. If the H20 queue is very long, you can force it onto the L20 instead by submitting as:
> ```bash
> sbatch -p l20-gpu jobs/fam129b_capz.job
> ```

### 7.4 What the job script actually does (for reference)

You do not need to edit these files. But here is what happens when you `sbatch` a job:

1. Slurm reserves one GPU and 10 CPU cores
2. The script creates input/output directories in `/work/NETID/`
3. It copies the JSON file to `/work/NETID/AF3_Input/`
4. It runs AlphaFold 3 inside a Singularity container using the Duke installation
5. Results are saved to `/work/NETID/AF3_Output/<job_name>/`
6. Log files go to `jobs/logs/` (one `.out` file and one `.err` file per job)

---

## 8. Monitoring your jobs

### 8.1 Check the queue

```bash
squeue -u NETID
```

Example output:

```
JOBID     PARTITION  NAME                 USER  ST  TIME  NODES
1234567   l20-gpu    AF3_fam129b_rab5a_g  abc12 R   0:15  1
1234568   l20-gpu    AF3_fam129b_rab5a_g  abc12 PD  0:00  1
```

Status codes:
- `PD` = pending (waiting for a GPU to be free)
- `R` = running
- `CG` = completing
- Nothing listed = the job finished (or failed)

### 8.2 Check if a job finished successfully

```bash
cat jobs/logs/slurm-AF3_fam129b_rab5a_gtp_1234567.out
```

The last two lines of a successful job look like:

```
The analysis is done.
Mon Sep 29 14:23:01 EDT 2026
```

### 8.3 Check for errors

```bash
cat jobs/logs/slurm-AF3_fam129b_rab5a_gtp_1234567.err
```

If this file is empty or contains only minor warnings, the job succeeded. Real errors usually say things like `ERROR`, `FileNotFoundError`, or `CUDA out of memory`.

### 8.4 Check if output files were created

```bash
ls /work/NETID/AF3_Output/fam129b_rab5a_gtp/
```

A successful run produces these files:

```
fam129b_rab5a_gtp_model.cif          ← the 3D structure
fam129b_rab5a_gtp_summary_confidences.json   ← the confidence scores
fam129b_rab5a_gtp_confidences.json
ranking_scores.csv
seed-1_sample-0/   seed-1_sample-1/   ... seed-1_sample-4/
```

---

## 9. Copying results back to your laptop

Once jobs finish, run this **on your laptop** (not on the cluster). Replace `NETID`:

```bash
rsync -av NETID@dkucc-login-01.rc.duke.edu:/work/NETID/AF3_Output/ /Users/jetyue04/af3/output/
```

This copies all output into the `output/` folder next to the pipeline code. The analysis scripts know to look there automatically.

---

## 10. Running the analysis pipeline

All analysis runs on your laptop. The pipeline has two parts: a quick scoring script, then a Jupyter notebook for visual review.

### 10.1 One-time environment setup

The first time only — installs all required software:

```bash
cd /Users/jetyue04/af3/af3-pipeline
conda env create -f environment.yml
```

This takes about 5–10 minutes. You only need to do this once.

To verify it worked:

```bash
conda activate af3-pipeline
python -c "import pandas, Bio, rdkit; print('OK')"
```

If it prints `OK`, you are ready.

### 10.2 Step 1 — Extract confidence scores from AF3 output

```bash
conda activate af3-pipeline
cd /Users/jetyue04/af3/af3-pipeline
python src/scoring/parse_af3.py
```

This reads the `summary_confidences.json` files from every completed AF3 run and writes a summary table to `results/af3_scores.csv`. It runs in a few seconds.

### 10.3 Step 2 — Open the validation notebook

```bash
conda activate af3-pipeline
cd /Users/jetyue04/af3/af3-pipeline
jupyter lab notebooks/pipeline_review.ipynb
```

A browser window will open automatically. The notebook has 6 sections:

| Section | What it does | Ready to run? |
|---------|-------------|---------------|
| **0 — Setup** | Loads libraries, checks paths | ✅ Always |
| **1 — Monomer QC** | Plots pLDDT confidence for each protein alone | ✅ With existing data |
| **2 — Control validation** | Checks that RAB5A-GTP scores better than RAB5A-GDP | ✅ With existing data |
| **3 — FAM129B complexes** | Scores and plots the novel protein pairs | After cluster jobs finish |
| **4 — Docking setup** | Extracts the binding interface, prepares the receptor | After Section 3 passes |
| **5 — Screening results** | Shows ranked small-molecule hits | After docking runs |

**How to run a section:** Click inside any code cell and press `Shift + Enter`. Run cells from top to bottom in order.

Each section ends with a printed verdict:
- `✓ PASS` — safe to continue
- `✗ FAIL` or `⚠ WARNING` — read the message carefully before going further

> **Do not proceed to FAM129B analysis if Section 2 fails.** The biology control must hold first.

### 10.4 Step 3 — Virtual screening (after Section 3 passes)

Once the notebook confirms the FAM129B complex structures are good, run the docking screen:

```bash
conda activate af3-pipeline
cd /Users/jetyue04/af3/af3-pipeline
python src/docking/screen.py
```

This takes 1–4 hours depending on the library size. It screens ~5,000 small molecules and writes `results/docking_hits.csv`.

### 10.5 Adding your own compounds to screen

Edit `ligands/custom.smi` — one compound per line, format `SMILES name`:

```
CC1=CC=C(C=C1)NC(=O)C  my_compound_001
c1ccc2c(c1)cccc2       my_compound_002
```

SMILES strings for most compounds can be found on PubChem or ChemDraw. Your compounds will automatically be included in the next screening run.

---

## 11. What the scores mean

### AlphaFold 3 scores

| Score | What it measures | Good threshold | Interpretation |
|-------|-----------------|----------------|---------------|
| `ptm` | Overall structure confidence (single protein) | > 0.7 | How reliable is the monomer structure |
| `iptm` | Interface confidence (protein complex) | > 0.75 = strong; 0.5–0.75 = weak | How confident AF3 is that these two proteins actually bind |
| `fraction_disordered` | % of residues with low confidence | < 0.3 ideally | Some disordered regions are expected; FAM129B naturally has a long disordered tail |
| `ranking_score` | Weighted combination; AF3's own ranking | Higher is better | Useful for comparing multiple runs of the same job |

**The critical check:** `iptm(GTP state) > iptm(GDP state)` for the same protein pair. This means AF3 is correctly detecting that the complex forms preferentially when RAB5A is in its active (GTP-bound) conformation — consistent with known biology.

### Docking scores

| Score | Meaning | Hit threshold |
|-------|---------|---------------|
| `vina_score` | AutoDock Vina binding energy (kcal/mol) | More negative = better binding |
| `vina_zscore` | Statistical significance vs. 20 random decoy compounds | ≤ −2.5 = likely hit |
| `gnina_score` | Deep-learning rescore (independent confirmation) | More negative = better |
| `consensus_flag` | `True` if compound scores well in BOTH Vina AND GNINA | `True` = high confidence hit |

**Top hits** are compounds with `consensus_flag = True` and the most negative `vina_zscore`. These are the candidates to order for experimental testing.

---

## 12. Full pipeline at a glance

```
YOUR LAPTOP                         DUKE CLUSTER
───────────────────────────────     ─────────────────────────────────────
 1. git pull                        
    (get latest code)               
                                    
 2. rsync laptop → cluster   ──────►  /work/NETID/af3-pipeline/
    (copy code to cluster)          
                                    
                                     3. ssh to cluster
                                        sbatch jobs/fam129b_rab5a_gtp.job
                                        sbatch jobs/fam129b_rab5a_gdp.job
                                        (... 9 more jobs)
                                    
                                     4. squeue -u NETID
                                        (monitor; each job ~1–3 hrs)
                                    
                                     5. Jobs complete →
                                        /work/NETID/AF3_Output/
                                    
 6. rsync cluster → laptop   ◄──────  rsync AF3_Output/ → local output/
    (copy results back)             
                                    
 7. conda activate af3-pipeline     
    python src/scoring/parse_af3.py 
    → results/af3_scores.csv        
                                    
 8. jupyter lab                     
    notebooks/pipeline_review.ipynb 
    Sections 0 → 1 → 2 → 3         
    (inspect, confirm controls pass)
                                    
 9. python src/docking/screen.py    
    → results/docking_hits.csv      
    (1–4 hrs on laptop)             
                                    
10. Notebook Section 5              
    Review top hits                 
    Shortlist 5–10 compounds        
    for wet-lab ordering            
```

---

## 13. Troubleshooting

### "ssh: Could not resolve hostname"
You are not on Duke VPN. Connect via Cisco AnyConnect first.

### "Permission denied (publickey)"
Your Duke password may have changed, or the SSH key needs to be set up. Contact Duke Research Computing at rc-support@duke.edu.

### Slurm job shows `PD` for a long time
The GPU queue is busy. `l20-gpu` and `h20-gpu` are shared cluster resources. Peak times (weekday afternoons) can have waits of several hours. Submit at night or on weekends for faster turnaround.

### "CUDA out of memory" in the `.err` log
The job ran out of GPU memory. Switch the job to the H20 partition (more VRAM):
```bash
sbatch -p h20-gpu jobs/your_job.job
```

### `python src/scoring/parse_af3.py` finds no output files
Either the jobs have not finished yet, or the rsync copy in Step 9 has not been run. Check with:
```bash
ls /Users/jetyue04/af3/output/
```
If the output folder is empty, run the rsync command from Section 9.

### Notebook cell shows a red error about missing module
The conda environment is not activated. In Terminal:
```bash
conda activate af3-pipeline
jupyter lab notebooks/pipeline_review.ipynb
```

### "version": 4 causes job failure
The Duke cluster's AlphaFold 3 installation only accepts `"version": 2` in the JSON. All existing JSONs are already set correctly. If you create a new JSON from an online template, change `"version": 4` to `"version": 2`.

---

*For questions about the code or pipeline logic, contact your graduate student. For cluster access, accounts, or quota issues, contact Duke Research Computing at rc-support@duke.edu.*
