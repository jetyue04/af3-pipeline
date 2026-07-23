# FAM129B / RAB5 / RAB7 / CapZ — AlphaFold 3 Screening

## Right now: manual walkthrough, not the agent

Before building any orchestration/agent (see the design docs in
[claude/](claude/) for that later phase), the immediate goal is to run the
AlphaFold 3 pipeline **by hand, one step at a time**, to actually understand
what each stage produces and how the SLURM job behaves. No automation yet —
just enough manual reps to trust the pipeline before wrapping code around it.

### Steps to do manually

1. **Fetch sequences from UniProt** (by hand, one `curl`/browser download at a
   time — no `fetch.py` yet):
   - RAB5A — P20339
   - RAB7A — P51149
   - FAM129B — Q96TA1
   - CAPZA1 — P13127
   - CAPZB — accession still TBD, look up the reviewed Swiss-Prot entry
2. **Hand-build one AF3 input JSON** for a single test pair (start with a
   known Tier 1 control, e.g. RAB5A + EEA1, not a FAM129B hypothesis — known
   biology first, so a bad result means "I set up the job wrong," not
   "the hypothesis is wrong").
3. **Submit it on SLURM** using [jobs/af3-test.job](jobs/af3-test.job) (or the
   `af3-test1.job` variant) and watch it run start to finish.
4. **Look at the raw output directory** — get familiar with what AF3 actually
   writes out (structure files, confidence metrics) before writing any
   parsing code for it.
5. **Manually read the scores** — pLDDT, PAE, ipTM — for that one job, and
   sanity-check them against what's expected for a real interaction.
6. Repeat steps 2–5 for a couple more Tier 1 controls. Only move on to a real
   FAM129B pair once the manual process feels understood end to end.

## Repo layout

- `jobs/` — SLURM `.job` scripts for submitting AF3 runs.
- `claude/` — planning docs (not active work): `pipeline.md` is the
  project-specific protocol/repo-design notes, `guide.md` is a more generic
  agent/docking blueprint. Both describe the *eventual* automated pipeline —
  reference them later, not now.
- `protein/` — sequence/structure data as it's fetched (currently empty,
  filling in as step 1 above happens).

## Later — not now

Once a handful of manual runs make sense, revisit `claude/pipeline.md`'s
proposed `src/` domain layout (`protein/`, `complex/`, `af3/`, `scoring/`,
`pipeline/`) and `claude/guide.md`'s agent-loop blueprint to start automating
the parts that turned out to be repetitive.
