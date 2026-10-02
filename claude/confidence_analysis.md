# AF3 Confidence Score Analysis & Hypothesis Status

**Date:** 2026-09-24

---

## All Jobs — Confidence Scores

| Job | iptm | ptm | pair iptm (A↔B) | ranking | disordered |
|-----|------|-----|-----------------|---------|------------|
| fam129b_niban2_monomer | — | 0.75 | — | 0.85 | 21% |
| rab5a_full_monomer | — | 0.76 | — | 0.86 | 21% |
| rabep1_rabaptin5_monomer | — | 0.23 | — | 0.53 | 60% |
| rab5a_full_gtp_mg | 0.97 | 0.80 | 0.97 | 1.04 | 21% |
| rab5a_full_gdp_mg | 0.97 | 0.77 | 0.97 | 1.03 | 21% |
| rab5a_gtp_mg_rabaptin5_complex | 0.70 | 0.34 | 0.68 | 0.82 | 40% |
| rab5a_gdp_mg_rabaptin5_complex | 0.69 | 0.32 | 0.68 | 0.84 | 43% |
| fam129b_rab5a_gtp_mg_complex | 0.28 | 0.60 | 0.16 | 0.44 | 20% |
| fam129b_rab5a_gdp_mg_complex | 0.33 | 0.64 | 0.28 | 0.49 | 20% |
| fam129b_capz_heterodimer_complex | 0.36 | 0.50 | 0.18 | 0.47 | 17% |

---

## Key Benchmarks

- **pair iptm > 0.68** = confident direct binder (Rabaptin5 benchmark)
- **pair iptm 0.16–0.28** = what FAM129B scores — 2–4x below benchmark
- **GTP iptm ≥ GDP iptm** = expected for a true GTP-state effector
- FAM129B shows **GDP > GTP** (0.28 vs 0.16) — opposite of effector pattern

---

## Interface Residues (5Å cutoff, GTP complex)

**RAB5A (chain B):**
- P-loop: Q20, K22
- Switch II: F57, L58, T59, K70, E72, W74
- Other: E50, I53, G54, A55, A56, R81, Y82, L85, M88, Y89

**FAM129B (chain A) — two clusters:**
- Cluster 1 (residues 359–378): T359, R362, D363, F366, K367, T370, D371, N373, L374, I377, N378
- Cluster 2 (residues 463–478): E463, C466, K467, Q470, E474, L477, K478

**Hotspot residues on RAB5A:** W74, F57 (hydrophobic Switch II anchor) + K22 (P-loop, near phosphate)
**Charged FAM129B residues pointing at phosphate:** R362, K367, K467, K478

---

## Hypothesis Status

### Original hypothesis: FAM129B is a RAB5A effector
- **Weakly supported** — pair iptm (0.16 GTP, 0.28 GDP) far below Rabaptin5 benchmark (0.68)
- GDP > GTP is the wrong direction for a GTP-state effector

### Alternative hypothesis: FAM129B is a RAB5A GAP
- **Consistent with data:**
  - GDP complex scores higher (product state more stable after hydrolysis)
  - Contacts land on P-loop (K22) and Switch II — exactly where GAP arginine finger inserts
  - FAM129B R362 is an arginine pointing toward the phosphate — candidate arginine finger
  - Low iptm fits transient catalytic interaction (GAPs bind briefly then release)
- **To verify:**
  1. Check FAM129B for TBC domain (InterPro/Pfam) — most RAB GAPs have one; non-TBC GAPs exist (e.g. RN-tre for RAB5)
  2. GTPase activity assay — does FAM129B accelerate RAB5A hydrolysis?
  3. Literature: does FAM129B overexpression reduce RAB5A-GTP levels or shrink early endosomes?

### Other possibilities
- Indirect interaction (requires scaffold, membrane, or co-factor)
- Weak/transient effector (Kd > 10 µM — AF3 struggles with these)
- Wrong domain — only a specific FAM129B domain may be responsible

---

## Visualizations Generated

- `results/fam129b_rab5a_gtp_complex.html` — interactive 3D, GTP state
- `results/fam129b_rab5a_gdp_complex.html` — interactive 3D, GDP state
- Both have 4Å / 5Å / 8Å distance filter buttons and per-residue tables

---

## Next Steps (pending PI discussion)

1. **Discuss GAP hypothesis with PI** before continuing to docking — drug strategy differs
2. Check FAM129B domain annotation (TBC domain?)
3. Wait for RAB7A complex results to see if pattern repeats
4. If proceeding to docking: best candidate receptor is RAB5A-GTP, chain B, Switch II patch
5. Run Stage 2 full analysis (BSA + PRODIGY ΔG) to get quantitative interface metrics
