# AF3 input templates

Standard AlphaFold 3 input JSON, in the shapes this project actually needs.
Copy the one that matches what you're testing, rename it, and fill in the
`REPLACE_WITH_*` placeholders — save it into [../inputs](../inputs), then
copy [../../jobs/template.job](../../jobs/template.job) and point it at
that file (see [jobs/README.md](../../jobs/README.md)).

- **`monomer.json`** — single protein, e.g. the QC-each-protein-alone step
  before any pairing (pipeline.md step 3).
- **`binary_complex.json`** — two proteins in one job, e.g. RAB5A + FAM129B
  (step 5).
- **`nucleotide_state.json`** — a RAB protein plus GTP + Mg²⁺ as ligands, for
  the GTP-bound ("on") state (step 4). Swap the ccdCode to `"GDP"` for the
  GDP-bound ("off") negative control.
- **`heterodimer_capz.json`** — CAPZA1 + CAPZB together, since CapZ is never
  modeled as one subunit alone (step 6).

## Field notes

- `id` — chain ID(s) for that entity. A list (e.g. `["A", "B"]`) makes
  multiple identical copies of the same sequence; each entity in `sequences`
  needs its own unique id(s).
- `sequences[].protein` / `.rna` / `.dna` / `.ligand` — one of these per
  entry. Ligands take either `ccdCodes` (PDB Chemical Component Dictionary
  codes, e.g. `GTP`, `GDP`, `MG`, `ATP`) or a `smiles` string, not both.
  `sequence` is required for protein/rna/dna.
  See "modifications" in the specs to include PTMs, if needed.
- `modelSeeds` — list of ints; AF3 runs the full pipeline once per seed.
  More seeds = more samples to compare, at the cost of runtime.
- `dialect` / `version` — leave as `"alphafold3"` / `2`. This is the input
  *schema* version, separate from the container version (the SLURM job uses
  `AF3_v3.0.1.sif`) — empirically, this cluster's AF3 build rejects anything
  above `2`, even though the public docs describe versions up to 4. That
  means no `description` field either (a v4-only addition) — don't add one
  unless the cluster's AF3 build gets upgraded and this is re-verified.
- Optional and not included here: `bondedAtomPairs` (explicit covalent bonds,
  e.g. to a modified residue), and precomputed `unpairedMsa`/`pairedMsa`/
  `templates` (skip these — the SLURM job already points at a full
  `AF3_databases` dir and runs MSA search itself).
