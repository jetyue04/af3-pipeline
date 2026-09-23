# Endosomal Trafficking Pipeline Blueprint: Dual Screening & Commercial IP Strategy

This document outlines the high-performance, legally compliant computational pipeline designed for screening small-molecule inhibitors targeting endosomal trafficking proteins (e.g., PD-L1, Rab5/Rab7 pathways). 

## 1. The Core Strategy: Exploiting the Legal Loophole
To protect the Principal Investigator's (PI) biotech startup from intellectual property contamination while maintaining high academic prestige, the lab will execute a **two-tier modeling strategy**:

*   **Phase 1 (Discovery Engine):** Driven entirely by **Boltz-2** (or Chai-1). These open-source tools match AlphaFold 3's accuracy but use highly permissive **MIT/Apache 2.0 licenses**. All startup patents and commercial assets remain 100% legally clean.
*   **Phase 2 (Academic Flex):** The official **AlphaFold 3 (AF3)** infrastructure is utilized *only at the very end* on the single winning compound. This generates the precise high-resolution figures required by top-tier journal reviewers (e.g., *Nature Cell Biology*, *BMC Biology*) within legal academic boundaries.

---

## 2. Approach A: The AI-First Pipeline (For Unknown / Mutated Proteins)
Use this workflow when the target endosomal protein has a rare sequence, heavy mutations, or has never been experimentally solved in the public database.

Use code with caution.[ Step 1: AI Co-Folding ] ────────> [ Step 2: AutoDock Vina ] ──────> [ Step 3: Wet-Lab Assay ]Feed target sequence +              Split the .cif file into          Order the physicaldrug candidate into Boltz-2.        receptor & ligand. Run            compound and test itGenerates flexible 3D shape.        Vina to check for clashes.        in live cancer cells.
### Step 1: Flexible AI Co-Folding & Native Affinity (Boltz-2)
*   **Action:** Feed the target protein sequence and the library of candidate drug SMILES strings directly into Boltz-2 on the cluster's GPUs.
*   **Purpose:** Builds the 3D protein structure from scratch and allows the binding pocket to dynamically "morph" around the drug atom-by-atom (induced fit). Boltz-2 predicts the structural coordinates and provides a native AI binding affinity score.

### Step 2: Hard Physical Audit (AutoDock Vina Evaluation)
*   **Action:** Split the Boltz-2 `.cif` output into separate protein and ligand files using a Python parser (Gemmi/Meeko). Pass them into AutoDock Vina with a search depth configuration of zero (`--exhaustiveness 0`).
*   **Purpose:** Vina acts as a strict physics engine to audit the AI's geometry. It checks electrostatic potentials at those exact coordinates and rejects any "hallucinated" poses where atoms physically overlap (steric clashes).

---

## 3. Approach B: The Physics-First Pipeline (For Known / Solved Proteins)
Use this workflow when you have access to an exact, experimentally solved 3D structure of the target protein (from the public Protein Data Bank) and want to screen a massive chemical library efficiently.

[ Step 1: CPU Screening ] ──> [ Step 2: GPU Co-Folding ] ──> [ Step 3: Hard Physics Check ]10,000 drugs run via          Top 500 run via Boltz-2       AutoDock Vina evaluatesAutoDock Vina (Rigid)         (Flexible AI Affinity)        final 3D atomic coordinates
### Step 1: High-Throughput CPU Screening (AutoDock Vina)
*   **Action:** Take the known PDB structure and screen your entire library (e.g., 10,000+ compounds) using traditional AutoDock Vina on the cluster's CPUs.
*   **Purpose:** Rapidly eliminates thousands of non-binders in seconds, filtering out obvious misses and drastically reducing the GPU compute requirements.

### Step 2: AI Structural Generation (Boltz-2)
*   **Action:** Extract the top 500 promising compounds that passed Step 1 and pass them to Boltz-2 on your GPU node.
*   **Purpose:** Generates a relaxed, flexible co-folded model of the complex, resolving the "induced fit" mechanics that rigid CPU screening ignores.

### Step 3: Final Physics Score (AutoDock Vina Evaluation)
*   **Action:** Split the resulting Boltz-2 `.cif` file and execute an AutoDock Vina score-only audit (`--exhaustiveness 0`).
*   **Purpose:** Computes the verified physical binding energy of the flexible, AI-optimized binding coordinates.

---

## 4. Definitive Biological Validation (Bypassing MD)
*   **Do we need MD Simulations (e.g., GROMACS/FEP)?** **No.** 
*   **Justification:** While physics-based simulations like Free Energy Perturbation (FEP) are a computational gold standard, they require massive hardware allocations. Because the lab is transitioning directly into **wet-lab assays**, *in vitro* cell culture data and *in vivo* animal models provide definitive verification. 
*   **The Workflow Action:** Use the consensus rank from your selected pipeline approach to pick the top 5–10 high-confidence drug leads, order the physical compounds, and move directly to live testing on cell cultures.

---

## 5. Mandatory Steps for Legal IP Protection & Auditing
To ensure a verifiable data paper trail for future venture capital or pharmaceutical due diligence, the pipeline must automate the following documentation:

1.  **Git Code Commit Logs:** Store the entire screening pipeline execution logic in a private repository. Timestamped code changes prove that `boltz predict` executed the discovery, not AF3.
2.  **Immutable System Manifests:** Configure Slurm batch scripts to print an execution receipt (`run_manifest.txt`) storing the exact system timestamp, cluster hostname, and the mathematical MD5 checksum of the Boltz container image.
3.  **Embedded Structural Metadata Tags:** Program the Python file parser to inject `REMARK 999` text headers directly into the top of every generated `receptor.pdbqt` and `ligand.pdbqt` file, explicitly stating the structure was forged under an open-source MIT license.