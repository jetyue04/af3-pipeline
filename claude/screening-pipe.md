# AI-Driven Virtual Screening Pipeline: Targeting the RAB5–FAM129B Interface

This document defines a hybrid Machine Learning (ML) and physics-based structural bioinformatics pipeline. The objective is to identify small-molecule competitive inhibitors that sit directly at the contact interface between **RAB5** and **FAM129B** to disrupt their interaction.

---

## 🏗️ Pipeline Architecture Overview

The workflow uses a four-phase funnel strategy to maximize geometric sampling accuracy while maintaining high screening throughput.

Use code with caution.[Phase 1: Boltz-2 Ensemble] ──> Stochastic 3D structural sampling of the PPI complex│[Phase 2: PPI Baseline]     ──> Calculate natural ΔG via PRODIGY / FoldX│[Phase 3: High-Throughput Vina]─> Rigid receptor docking of targeted small-molecule library│[Phase 4: Triage & Rescore] ──> Convert scores to Z-scores; feed top hits back to Boltz-2
---

## 🛠️ Step-by-Step Implementation

### Phase 1: Protein-Protein Co-folding & Ensemble Generation
Because protein interfaces are dynamic ("breathing"), a single static structure is insufficient. We use **Boltz-2** to sample the conformational space.

1. **Stochastic Sampling**: Run Boltz-2 with the amino acid sequences of both RAB5 and FAM129B using multiple random seeds (`--seed`) to generate a pool of 40 structural outputs.
2. **Quality Triage**: Parse the output metadata JSON files. Retain only models with an interface Predicted Template Modeling (**ipTM**) score higher than `0.5`.
3. **Structural Superimposition**: Select the highest-scoring model as a reference. Use Biopython (`Bio.PDB`) or `MMalign` to structurally align the remaining ensemble structures based on the rigid backbone of the **RAB5** target protein.

### Phase 2: Mapping the Target Interface & Baseline Calculation
We map the physical space where the two proteins touch and calculate the natural binding energy that the small molecules must compete against.

1. **Calculate Natural Δ G**: Feed your aligned Boltz-2 ensemble structures into **PRODIGY** or **FoldX**. Record the average natural interface free energy (\(\Delta G_{\text{natural}}\)) in kcal/mol. This serves as your therapeutic baseline.
2. **Isolate the Receptor**: Strip the FAM129B protein chain out of the PDB files, leaving an isolated **RAB5** structure with its interaction face completely exposed.
3. **Define the Search Space**: Measure the average X, Y, and Z coordinates of the contact residues (focusing near the critical *FAM129B K83 site*). Use these coordinates to center your AutoDock Vina grid box, applying a buffered dimensions layout (e.g., 24 × 24 × 24 Å).

### Phase 3: High-Throughput Virtual Screening via AutoDock Vina
We use **AutoDock Vina** strictly for **normal (rigid receptor) docking**. This maintains optimal screening throughput while the incoming small molecules retain full flexibility.

1. **Library Preparation**: Curate a **PPI-focused chemical library** (such as Enamine or Life Chemicals PPI sets). Do not use general chemical libraries, as PPI interfaces require larger, flatter, aromatic structures. Convert the library from SMILES to `PDBQT` using RDKit and Meeko.
2. **Configure Vina**: Write a centralized `config.txt` targeting the isolated RAB5 contact interface:
   ```text
   receptor = rab5_isolated.pdbqt
   center_x = [Calculated_Interface_X]
   center_y = [Calculated_Interface_Y]
   center_z = [Calculated_Interface_Z]
   size_x = 24
   size_y = 24
   size_z = 24
   exhaustiveness = 8
   ```
3. **Execution Loop**: Run the automated screening loop across your ligand dataset:
   ```bash
   vina --config config.txt --ligand ligand_i.pdbqt --out out_i.pdbqt --log log_i.txt
   ```

### Phase 4: Z-Score Triage and Boltz-2 Rescoring
To prevent mathematical comparison errors between different software engines, normalization must be applied.

1. **Z-Score Normalization**: Run a small control set of 500 random decoy molecules to calculate the scoring average and standard deviation of your specific receptor model setup. Convert all raw Vina outputs to statistical Z-scores:
   $$\text{Z-Score} = \frac{\text{Ligand Vina Score} - \text{Decoy Average}}{\text{Standard Deviation}}$$
2. **Hit Filtering**: Isolate any compounds achieving a high negative Z-score (e.g., $Z \le -3.0$).
3. **The "Sandwich" Validation**: Take the top 100 3D ligand output coordinates (`out_i.pdbqt` converted back to PDB/SDF) and feed them back into the **Boltz-2 Binder Confidence Module**. If a ligand yields both a high negative Vina Z-score and a high Boltz-2 confidence classification score, it is selected for wet-lab ordering.

---

## 🧪 Wet-Lab Experimental Validation

The final computational hit list (top 10–20 consensus compounds) must be validated using low-throughput physical assays:

1. **Direct Binding (SPR / TSA)**: Purify raw recombinant RAB5. Use Surface Plasmon Resonance (SPR) or a Thermal Shift Assay (TSA) to prove the small molecules physically bind to the protein in isolation.
2. **Interface Disruption (AlphaScreen / TR-FRET)**: Tag RAB5 and FAM129B with donor/acceptor fluorescent beads. Introduce the drug candidates. A successful blocker will physically separate the proteins, causing a measurable drop in fluorescent signal.
3. **In Vivo Cellular Phenotype**: Treat living cells with the compound. Perform confocal microscopy tracking to observe if early endosome vesicle fusion (governed by RAB5) shifts or normalizes, confirming functional pathway disruption.