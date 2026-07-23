# Filled AF3 input JSONs

Ready-to-run AlphaFold 3 inputs, built from the templates in
[../templates](../templates) and the sequences in [../sequences](../sequences).
Copy whichever one to `AF3_Input/alphafold_input.json` on the cluster to run
it via [jobs/](../../jobs/).

- `rab5a_monomer.json` — plain full-length RAB5A (P20339), no ligand — the
  apo/nucleotide-free baseline (step 3).
- `rab5a_full_gtp.json` — full-length RAB5A (P20339) + GTP + Mg²⁺, the "on"
  state (pipeline.md step 4).
- `rab5a_full_gdp.json` — same, with GDP instead of GTP — the "off" state,
  used as a negative control (real effectors should prefer GTP).
- `fam129b_monomer.json` — FAM129B (Q96TA1) monomer, for the standalone QC
  pass before any pairing (step 3).
- `rabaptin5_monomer.json` — Rabaptin-5 / RABEP1 (Q15276) monomer, a known
  RAB5 effector used as a Tier 1 positive-control partner.
- `rab5a_rabaptin5_gtp.json` — binary complex: full-length RAB5A + Rabaptin-5
  + GTP + Mg²⁺. The Tier 1 positive-control pairing (step 5) — Rabaptin-5 is
  a real effector, so it should dock onto the GTP ("on") state.
- `rab5a_rabaptin5_gdp.json` — same pairing with GDP instead of GTP, the
  built-in negative control (a real effector should bind this one worse).

All model seeds are `[1]` — bump `modelSeeds` to `[1, 2, ...]` for more
samples once you're past the single-seed sanity-check stage.

**Note:** UniProt's current gene symbol for FAM129B is `NIBAN2` (it was
renamed) — accession Q96TA1 is still the right one, just filed under the new
name.
