# SLURM jobs

- `af3-test.job` — the **template**. Copy it for any new protein/complex,
  fill in the two `REPLACE_WITH_*` spots (job name, input JSON filename),
  and it's ready to submit. Also documents the cluster's GPU partitions.
- `af3-test1.job` — earlier scratch variant, expects the input JSON to
  already be sitting at `/opt/apps/community/alphafold3/AF3_Input/alphafold_input.json`.
  Kept for reference; new jobs should start from `af3-test.job` instead.
- One filled-in job per input in [../protein/inputs](../protein/inputs) —
  each is `af3-test.job` with the placeholders substituted, and a header
  comment noting exactly what was filled in. Each copies its matching JSON
  from this repo into `$INPUT/alphafold_input.json` before running, so you
  don't have to hand-edit anything per submission:
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

Output lands in `/work/$USER/AF3_Output`, and the slurm stdout/stderr logs
(`slurm-AF3_*_<jobid>.out/.err`) get written wherever you ran `sbatch` from.
