# SLURM jobs

- `template.job` — the **template**. Copy it for any new protein/complex,
  fill in the two `REPLACE_WITH_*` spots (job name, input JSON filename),
  and it's ready to submit. Also documents the cluster's GPU partitions.
- One filled-in job per input in [../protein/inputs](../protein/inputs) —
  each is `template.job` with the placeholders substituted, and a header
  comment noting exactly what was filled in. Each copies its matching JSON
  from this repo into the shared input dir before running, so you don't
  have to hand-edit anything per submission, and **all of them can be
  `sbatch`'d at the same time** without clobbering each other:
  - `rab5a_monomer.job`
  - `rab5a_full_gtp.job`
  - `rab5a_full_gdp.job`
  - `fam129b_monomer.job`
  - `rabaptin5_monomer.job`
  - `rab5a_rabaptin5_gtp.job`
  - `rab5a_rabaptin5_gdp.job`

## GPU partitions

This cluster has NVIDIA L20 (`l20-gpu`) and H20 (`h20-gpu`) nodes. Confirm
exact partition names/current load with `sinfo` before submitting — public
specs for the two GPUs (not this cluster's exact node counts/queue times):

| | L20 | H20 |
|---|---|---|
| VRAM | 48GB GDDR6 | 96GB HBM3 |
| Memory bandwidth | ~864 GB/s | ~4 TB/s |
| Architecture | Ada Lovelace | Hopper (NVLink) |

**H20 is the better GPU overall** — 2x the VRAM and ~4-5x the memory
bandwidth of L20, which is what AF3's attention layers (Evoformer/
Pairformer) are bottlenecked on, since they scale with sequence length. If
`h20-gpu` isn't queued/contended, default to it.

- **Use `h20-gpu`** whenever it's available — especially for multi-chain
  complexes, longer sequences, or deep MSA/templates, e.g. the
  RAB5A+Rabaptin-5 pairs and anything with FAM129B (726 residues) once it's
  paired up. These are the jobs where the extra VRAM/bandwidth matters most.
- **Use `l20-gpu`** as the fallback when H20 nodes are busy, or for single
  monomers and quick sanity-check runs where the weaker GPU is still plenty
  — keeps H20 nodes free for jobs that actually need them.

Switch partitions without editing the file via `sbatch -p <partition> ...`
(overrides the `#SBATCH -p` line in the script). The monomer jobs default to
`l20-gpu`; `rab5a_rabaptin5_gtp.job` and `rab5a_rabaptin5_gdp.job` (2-chain
complexes) are set to `h20-gpu`.

## Usage

These assume the repo is cloned into the SLURM work directory (default path
`/work/$USER/af3-pipeline`). If it lives somewhere else, set `REPO_DIR`
before submitting:

```bash
sbatch jobs/rab5a_monomer.job
# or, if the repo isn't at /work/$USER/af3-pipeline:
REPO_DIR=/path/to/af3-pipeline sbatch jobs/rab5a_monomer.job
```

## Input/output layout

- Input: `/work/$USER/AF3_Input/<job_name>.json` — a single **shared**
  directory. Each job just copies its JSON in under its original filename
  from `protein/inputs/`; no per-job subfolder needed, because
  `run_alphafold.py` is called with `--json_path=/root/af_input/<job_name>.json`
  pointing at that one specific file — not `--input_dir`, which would
  process every file sitting in the directory. That's what makes it safe
  for multiple jobs to share one `AF3_Input` dir and run concurrently.
- Output: `/work/$USER/AF3_Output/<job_name>/` (own subfolder per job, via
  `JOB_NAME` in the script) — e.g. `/work/$USER/AF3_Output/rab5a_rabaptin5_gtp/`.

Slurm stdout/stderr logs (`slurm-AF3_*_<jobid>.out/.err`) get written
wherever you ran `sbatch` from, same as before.

## Running everything at once

```bash
cd /work/$USER/af3-pipeline
for j in jobs/*.job; do
  [[ "$j" == *template.job ]] && continue
  sbatch "$j"
done
```

Safe to fire off in one go: each job only ever touches its own JSON via
`--json_path`, and writes to its own `AF3_Output/<job_name>/` subfolder.
