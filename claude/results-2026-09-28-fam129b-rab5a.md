# FAM129B × RAB5A AF3 results — session notes (2026-09-28)

Status: **analyzed, inconclusive**. AF3 does not confidently support FAM129B–RAB5A
binding in either nucleotide state. Do NOT proceed to docking on these models.
Continued next session — open questions and next steps at bottom.

---

## What was run

Existing pipeline scripts only (`conda env: af3-pipeline`):

```bash
cd /Users/jetyue04/af3/af3-pipeline
AF3_OUTPUT_DIR=/Users/jetyue04/af3/output python src/scoring/parse_af3.py
```

- `src/scoring/parse_af3.py` → `results/af3_scores.csv`
- `src/scoring/interface.py` (`compute_interface`, `measure_switch_distances`) on
  chains A/B of each top-ranked `*_model.cif`
- `src/scoring/prodigy.py` (`run_prodigy`) on the same chains
- Per-sample check: all 5 `seed-1_sample-*/model.cif` per job (note: sample CIFs
  are named `model.cif`, NOT `*_model.cif`)

Chain mapping (both jobs): A = FAM129B (746 aa), B = RAB5A (215 aa),
C = GDP/GTP, D = Mg²⁺. Results table: `results/fam129b_binding_summary.csv`.

## Headline numbers

| Metric (FAM129B↔RAB5A) | GDP job | GTP job | Rabaptin5 control (known binder) |
|---|---|---|---|
| ipTM (complex) | 0.33 | 0.28 | 0.69 (GDP) / 0.70 (GTP) |
| **Pair ipTM (A↔B)** | **0.28** | **0.16** | — |
| **Pair PAE min (A↔B, Å)** | **19.9** | **24.9** | — |
| BSA (Å²) | 1604 | 1352 | 2055 (GDP) / 1695 (GTP) |
| Contact pairs ≤5 Å | 49 | 43 | 58 (GDP) / 46 (GTP) |
| PRODIGY ΔG (kcal/mol) | −9.37 | −8.86 | −12.61 (GDP) / −9.94 (GTP) |
| PRODIGY Kd | ~134 nM | ~315 nM | 0.56 nM – 51 nM |

Individual proteins fold fine in both jobs (chain pTM 0.73–0.75); RAB5A↔nucleotide
interfaces tight (0.93–0.95). Failure is specific to the protein–protein interface.

## Interpretation (agreed direction)

- Pair ipTM 0.16–0.28 is far below even the "weak" band (0.5–0.75); PAE 20–25 Å
  means AF3 has no idea of the relative orientation. Classic signature: both
  proteins folded correctly, docking did not.
- **GTP-vs-GDP verdict: backwards.** Real RAB5 effectors prefer GTP ("on" state).
  GDP scored ≥ GTP on every metric. No nucleotide-state specificity.
- PRODIGY ΔG (~−9 kcal/mol) is geometry-only — it assumes the modeled interface
  is real. Given the low ipTM/PAE, treat as "if the interface were real"
  estimate, NOT evidence of binding. Do not quote these ΔGs as affinity data.
- Per-sample consistency: GDP samples agree with each other (44–49 Å COM, 45–56
  contacts) = consistently weak; GTP samples disagree (27–60 contacts,
  1125–2009 Å² BSA) = no reproducible placement. "GDP > GTP" = GTP inconsistent,
  not GDP good.

## Competing explanations (not yet resolved)

**A. Prediction is bad / uninformative** — weak AF3 ≠ proof of non-binding:
1. Full-length FAM129B (20% disordered) — guide.md step 3 already anticipated
   fragmenting FAM129B once disorder is mapped. Highest-value fix.
2. Membrane context missing (guide.md step 11): RAB5 is prenylated; effector
   binding may be membrane-dependent.
3. AF3 known to miss weak/transient/fuzzy interfaces.

**B. Alternative biology (if taken at face value):**
- H1: FAM129B is not a *direct* RAB5 effector — interaction indirect (3rd
  protein), or direct partner is RAB7A or CapZ instead. Note:
  `fam129b_capz_heterodimer_complex` also weak (ipTM 0.36) — needs same scrutiny.
- H2: FAM129B binds GDP state = negative regulator (GDI-like). Weakly supported
  (GDP ≥ GTP) but within noise; speculative.
- H3: Real interface that AF3 cannot see (disorder/membrane/PTM-mediated).

## Next steps (priority order)

1. **Trim FAM129B construct** — pull monomer pLDDT (`fam129b_niban2_monomer`
   output), identify folded core, re-run co-folding vs RAB5A-GTP/GDP with trimmed
   construct (needs cluster/VPN + SLURM jobs in `jobs/`).
2. **Wet-lab tiebreaker** — nucleotide-state pull-down/co-IP (GTPγS-locked vs
   GDP-locked RAB5A). Definitive for effector-vs-not; AF3 was only triage.
3. Optional: cross-docking (HADDOCK/ClusPro, guide.md steps 8–9) — LOW value until
   a reproducible AF3 interface exists; skip for now.

## Loose ends for tomorrow

- [ ] **Notebook never executed + naming bug.** `notebooks/pipeline_review.ipynb`
      references jobs as `fam129b_rab5a_gtp` / `fam129b_rab5a_gdp` / `fam129b_capz`
      but actual dirs are `fam129b_rab5a_gtp_mg_complex`, `fam129b_rab5a_gdp_mg_complex`,
      `fam129b_capz_heterodimer_complex`. All cells unexecuted (no outputs).
      Fix names, then run Sections 0–3 (stop before docking Section 4).
      Section 3's own verdict for FAM129B×RAB5A would print "✗ WEAK / NO
      INTERACTION" (requires GTP > GDP on BSA+contacts; observed opposite).
- [ ] `rab7a_*` AF3 jobs (from README) — check if ever submitted; needed for H1.
- [ ] RAB7A/RAB5A G-domain constructs (guide.md step 1) — trimmed RAB constructs
      were also part of the plan; full-length RAB5A used here.
- [ ] Driver script from today was in `/tmp/run_interface_analysis.py` — recreate
      from this note if needed (imports `src.scoring`, chain IDs A/B hardcoded).
