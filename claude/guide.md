# FAM129B / RAB5 / RAB7 / CapZ — AlphaFold 3 Screening Pipeline

Context doc for Claude Code. Summarizes the research background, the protocol
we're implementing, and the repo/pipeline design decisions made so far.

## Background: what this research is

This project supports Jianbo Yue's lab. The scientific question: does **FAM129B**
physically interact with **RAB5**, **RAB7**, and/or **CapZ** (CAPZA1 + CAPZB)?

- **RAB5 / RAB7**: small GTPase "traffic control" proteins that label stages of
  endocytosis. RAB5 marks early endosomes (just internalized); RAB7 marks late
  endosomes (heading toward degradation/recycling). Both act as molecular
  switches: **GTP-bound = "on"** (active, recruits effector proteins),
  **GDP-bound = "off"**. Real effectors should preferentially bind the GTP
  ("on") state — this is used as a built-in control throughout the protocol.
- **FAM129B**: a less-characterized protein suspected of being a novel RAB5/RAB7
  *effector* (a helper protein recruited by the active, GTP-bound RAB).
- **CapZ**: an actin-filament-capping protein complex — always a heterodimer of
  **CAPZA1 + CAPZB**, never modeled as one subunit alone.

**Central hypothesis**: FAM129B may bridge RAB-mediated trafficking (RAB5/RAB7)
with cytoskeletal regulation (CapZ) — a potentially novel mechanistic link. This
has translational relevance to receptor trafficking pathways involving
**integrins** and **PD-L1**, both cancer-relevant.

**Why AlphaFold**: Traditionally, testing protein-protein interactions requires
slow, expensive wet-lab work (binding assays, crystallography/cryo-EM).
AlphaFold 3 predicts multi-protein complex structures computationally, acting as
a fast hypothesis-generation/triage step — flagging which pairs are worth
prioritizing for real experiments. **Critical caveat repeated throughout the
protocol: a confident AlphaFold prediction is a hypothesis, not proof.** It
cannot confirm an interaction happens in living cells, cannot model
biomolecular condensates, and cannot test drug/compound effects (e.g., "6J-1").

## Protocol summary (Jianbo Yue's document)

1. **Prepare protein sequences** — pull from UniProt: RAB5A (P20339), RAB7A
   (P51149), FAM129B (Q96TA1), CAPZA1 (P13127), CAPZB (accession TBD — look up
   reviewed/Swiss-Prot entry). For RAB5A/RAB7A, prepare both full-length and a
   trimmed **G-domain** construct (core folded domain, ~residues 1–185, minus
   the membrane-anchoring tail).
2. **Prioritized interaction matrix** — Tier 1: known RAB5/RAB7 partners
   (EEA1, RILP, etc.) used as **positive controls** to validate AlphaFold's
   reliability on established biology. Tier 2: the actual novel hypotheses
   (FAM129B × RAB5/RAB7/CapZ).
3. **Predict individual (monomer) structures first** — QC each protein alone
   (confidence scores, consistency across AF3's 5 samples) before pairing.
   FAM129B likely needs fragmenting once disordered vs. structured regions are
   identified.
4. **Model RAB5/RAB7 nucleotide states** — GTP and GDP separately, with Mg²⁺
   cofactor included. GTP ("on") state prioritized as the biologically relevant
   one.
5. **Predict binary complexes: RAB + FAM129B** — cross full-length/G-domain ×
   full/fragmented FAM129B × GTP/GDP. GDP comparison is a built-in negative
   control (real effectors should prefer GTP).
6. **Predict FAM129B + CapZ** — as full heterodimer, plus each CAPZ subunit
   alone (to help localize the binding interface, even though only the
   heterodimer is biologically realistic).
7. **Score/evaluate credibility** — record pLDDT, PAE, ipTM, contact counts,
   consistency across the 5 models. Verdict table translates score
   combinations into strong/weak/unstable calls. Note: AF3 can't detect
   biomolecular condensate-type interactions (different from fixed-structure
   binding).
8–9. **Independent docking cross-check** — re-test AF3-flagged interfaces with
   separate docking software (HADDOCK, ClusPro) since AF3's method isn't
   classic docking. Pair-specific interface checks (e.g., does FAM129B occupy
   CapZ's actin-capping surface, suggesting competition, or a different site?).
10. **Ternary/quaternary complexes** — only after binary pairs look solid,
    model 3–4 proteins together to distinguish possible mechanistic models
    (FAM129B as a bridge vs. CapZ binding RAB5 directly vs. competition for the
    same FAM129B site).
11. **Membrane context sanity check** — RAB proteins are normally
    membrane-anchored; isolated computational models can miss steric
    impossibilities that membrane geometry would rule out.
12. **Experimental validation** — the wet-lab phase: co-IP/pull-downs, binding
    assays, mutagenesis, cell imaging, functional assays (e.g., receptor
    recycling). This is what actually proves an interaction — AF3 only
    prioritizes candidates for this stage.

**Recommended first jobs**: a prioritized list of ~13 initial computational
jobs, starting with monomers and known controls before novel FAM129B pairs.

## Compute environment

- Running **AlphaFold 3 locally** on a **SLURM-based cluster** (this project),
  which requires **JSON input** (not the AlphaFold Server web UI, which takes
  FASTA/plain sequence paste — not applicable here).
- Separately has access to **NRP Nautilus** (Kubernetes-based) — prior
  AlphaFold experience there was via a ColabFold Docker image, with working
  knowledge of PVC setup, job YAML manifests, and GPU resource requests. Not
  the target environment for this project, but relevant prior experience.
- Also has access to a **SLURM cluster at Duke Kunshan University (DKU)** —
  separate system, not confirmed as the one used for this project.

## Repo structure decision

Chose a **domain-based** `src/` layout over a stage-numbered one, specifically
because the end goal is an **automated screening pipeline / agent**, not a
manually-run one-off protocol.

**Rationale:**
- An agent/orchestrator needs composable, independently callable functions
  (`fetch_sequence()`, `build_complex_json()`, `score_prediction()`) with real
  return values to branch on — not monolithic stage scripts that weld pipeline
  order to biology logic.
- Stage-numbered scripts assume a fixed linear order. A real screening pipeline
  needs to re-score old predictions under new criteria, add proteins without
  regenerating everything downstream, or selectively retry ambiguous cases —
  domain folders decouple "what kind of operation" from "what order it was
  first run in."

```
yue-fam129b-af3/
├── README.md
├── environment.yml
├── config/
│   ├── proteins.yaml            # UniProt IDs, construct boundaries (full vs G-domain)
│   ├── interaction_matrix.yaml  # Tier 1 (controls) vs Tier 2 (hypotheses) pairs
│   └── af3_settings.yaml        # model seeds, num samples, MSA settings
│
├── data/
│   ├── raw_sequences/           # untouched FASTA from UniProt, one per accession
│   ├── constructs/              # derived FASTA after trimming (e.g. G-domain)
│   └── ligands/                 # GTP/GDP + Mg2+ definitions for nucleotide states
│
├── src/
│   ├── __init__.py
│   ├── protein/
│   │   ├── fetch.py             # UniProt REST pulls
│   │   ├── construct.py         # trimming/fragmenting sequences
│   │   └── nucleotide.py        # GTP/GDP + Mg2+ state handling
│   │
│   ├── complex/
│   │   ├── pairs.py             # binary complex job generation (protocol steps 5-6)
│   │   └── ternary.py           # higher-order complex generation (step 10)
│   │
│   ├── af3/
│   │   ├── job_builder.py       # config -> AF3 input JSON
│   │   └── submit.py            # SLURM sbatch submission (or local run)
│   │
│   ├── scoring/
│   │   ├── parse.py             # extract pLDDT/PAE/ipTM from AF3 output
│   │   └── verdict.py           # step 7 rulebook -> strong/weak/unstable labels
│   │
│   ├── docking/
│   │   └── cross_validate.py    # HADDOCK/ClusPro cross-checks (steps 8-9)
│   │
│   └── pipeline/
│       ├── orchestrator.py      # DAG/state logic: what needs to run given current results
│       ├── state.py             # manifest of what's fetched/built/predicted/scored
│       └── cli.py               # entrypoint, e.g. `python -m pipeline run --stage scoring`
│
├── slurm/
│   └── job_templates/           # .sbatch templates, parameterized
│
├── inputs/
│   └── af3_json/                 # generated AF3 input JSONs
├── outputs/                      # gitignored — raw AF3 predictions are large/regenerable
├── results/
│   ├── scores_summary.csv        # tracked — the actual scientific record
│   └── figures/
├── notebooks/
└── docs/
    └── protocol.md                # Jianbo Yue's original protocol, versioned for reference
```

**Key points:**
- `outputs/` is gitignored (large, regenerable); `results/` (parsed scores) is
  tracked in git since that's the actual evolving scientific record.
- `state.py` is the piece that matters most for the "eventually an agent"
  goal — a manifest (JSON or SQLite) tracking what's been fetched/built/
  predicted/scored, with what config, so an orchestrator can decide what to
  run next instead of blind re-runs or silent skips.
- CAPZB's UniProt accession still needs to be looked up (reviewed/Swiss-Prot
  entry) and added to `config/proteins.yaml`.

## Open next steps (not yet started)
- Fetch and verify all sequences from UniProt REST API (`fetch.py`).
- Draft `config/proteins.yaml` and `config/interaction_matrix.yaml` schemas.
- Build `af3/job_builder.py` (config → AF3 JSON) and a SLURM `.sbatch` template.
- Decide on state-tracking format for `pipeline/state.py`.