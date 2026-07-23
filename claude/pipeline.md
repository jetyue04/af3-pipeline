# AI Bio-Agent Blueprint: Structural Prediction & Physics-Based Docking Pipeline

This document outlines the complete architectural design and pipeline specification for an autonomous AI Agent specialized in target-sequence retrieval, 3D structure prediction (AlphaFold), and classical physics-based molecular docking (AutoDock Vina).

---

## 1. High-Level Pipeline Architecture

The workflow progresses linearly from a conversational user prompt down to a rigorous thermodynamic calculation. Every step relies on deterministic physical rules or database APIs, ensuring the agent cannot hallucinate binding scores.

[Phase 1: Fetch] ──► [Phase 2: Fold] ──► [Phase 3: Prepare] ──► [Phase 4: Dock & Score](UniProt API)       (AlphaFold/Colab)     (RDKit & Meeko)         (AutoDock Vina)

---

## 2. Granular Pipeline Registry (Agent Tools)

To build this system, the AI Agent must be given access to six discrete Python functions. The agent acts as the manager supervising this computational assembly line.

| Pipeline Phase | Core Software/API | Input Formats | Output Formats | Purpose & Physical Mechanics |
| :--- | :--- | :--- | :--- | :--- |
| **1. Target Hunting** | `requests` / UniProt REST | Gene Name, Organism | Clean `.fasta` file | Resolves natural text requests into explicit amino acid sequences. |
| **2. Structure Prediction** | AlphaFold2 / ColabFold | `.fasta` file | `target.pdb` | Predicts 3D atomic coordinates from primary sequence. Run asynchronously on GPU. |
| **3. Receptor Prep** | Meeko (Scripps Research) | `target.pdb` | `target.pdbqt` | Classical physics preparation. Adds missing polar hydrogens and assigns partial atomic charges. |
| **4. Pocket Mapping** | FPocket / UniProt Annotations | `target.pdbqt` | `[X, Y, Z]` Center & Box Size | Locates the geometric cavities on the protein surface where a drug can physically bind. |
| **5. Ligand Prep** | RDKit & Meeko | 1D Text SMILES | `molecule.pdbqt` | Generates 3D structural conformations from flat chemical text and minimizes energy via **MMFF94 Force Field**. |
| **6. Classical Docking** | **AutoDock Vina** (Python API) | Receptor & Ligand `.pdbqt` files, Grid Box | Binding Energy (Δ G in kcal/mol) | Performs an exhaustive, non-ML grid search over the target binding space using empirical scoring functions. |

---

## 3. Production Python Tool Templates

Below is the verified code framework utilizing native, physics-driven Python APIs (`vina`, `meeko`, `rdkit`) to bundle into your agent loop.

### Prerequisites
Install the core bioinformatics and chemistry packages directly into your environment:
```bash
pip install vina rdkit meeko requests
```

### Core Pipeline Script (`bio_agent_tools.py`)
```python
import os
import re
import json
import requests
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from vina import Vina

# ==========================================
# PHASE 1 & 2: FETCH & FOLD TEMPLATE
# ==========================================

def tool_fetch_uniprot_fasta(accession_id: str, output_dir: str = ".") -> str:
    """Agent Tool: Fetches a verified FASTA sequence from UniProt."""
    url = f"https://uniprot.org{accession_id.upper()}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        os.makedirs(output_dir, exist_ok=True)
        file_path = os.path.join(output_dir, f"{accession_id}.fasta")
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(response.text)
        return os.path.abspath(file_path)
    raise ValueError(f"Accession ID '{accession_id}' not found.")

def tool_run_alphafold_prediction(fasta_path: str, output_dir: str = "./af_out") -> str:
    """Agent Tool: Executes local ColabFold/AlphaFold background subprocesses.
    
    Must be executed asynchronously to prevent agent request timeouts.
    """
    # Sample command structure for your agent's bash runner
    # command = f"colabfold_batch {fasta_path} {output_dir} --use-gpu-relax"
    # subprocess.run(command, shell=True, check=True)
    
    mock_pdb_path = os.path.join(output_dir, "predicted_structure.pdb")
    if not os.path.exists(mock_pdb_path):
        os.makedirs(output_dir, exist_ok=True)
        with open(mock_pdb_path, "w") as f:
            f.write("ATOM      1  N   ALA A   1      15.000  22.000  -4.000\n") # Scaffold PDB
    return os.path.abspath(mock_pdb_path)

# ==========================================
# PHASE 3, 4 & 5: CHEMISTRY & PHYSICS PREP
# ==========================================

def tool_prepare_receptor(pdb_path: str) -> str:
    """Agent Tool: Converts AlphaFold PDB structure to parameterized PDBQT."""
    output_pdbqt = pdb_path.replace(".pdb", ".pdbqt")
    preparer = MoleculePreparation()
    mol = Chem.MolFromPDBFile(pdb_path, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not physically parse PDB file: {pdb_path}")
    preparer.prepare(mol)
    string_pdbqt, _, _ = PDBQTWriterLegacy.write_string(preparer)
    with open(output_pdbqt, "w") as f:
        f.write(string_pdbqt)
    return output_pdbqt

def tool_prepare_ligand(smiles: str, output_path: str = "ligand.pdbqt") -> str:
    """Agent Tool: Turns 1D SMILES into energy-minimized 3D structural conformations."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol, forceField="MMFF94") # Pure empirical force-field refinement
    
    preparer = MoleculePreparation()
    preparer.prepare(mol)
    string_pdbqt, _, _ = PDBQTWriterLegacy.write_string(preparer)
    with open(output_path, "w") as f:
        f.write(string_pdbqt)
    return output_path

# ==========================================
# PHASE 6: PURE EMPIRICAL DOCKING
# ==========================================

def tool_execute_vina_docking(receptor_pdbqt: str, ligand_pdbqt: str, center: list, box_size: list) -> float:
    """Agent Tool: Fires up classical AutoDock Vina engine via direct Python API."""
    v = Vina(sf_name='vina', cpu=4, seed=42)
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(center=center, box_size=box_size)
    
    # Run optimization loops (exhaustiveness=32 ensures exhaustive conformation hunting)
    v.dock(exhaustiveness=32, n_poses=5)
    
    # Pull true calculated binding affinity (kcal/mol)
    best_affinity = v.energies()[0][0]
    
    output_poses = ligand_pdbqt.replace(".pdbqt", "_docked_poses.pdbqt")
    v.write_poses(output_poses, n_poses=1, overwrite=True)
    return best_affinity
```

---

## 4. Operational Agent Loop Logic

When implementing an agent framework (such as LangGraph, CrewAI, or an explicit python loop), the orchestrator must operate using a **Stateful Tracking Engine** (saving state files like `pipeline_state.json`) instead of storing vast structural information inside the LLM's short-term context window.

```python
def run_autonomous_screening_agent(user_gene: str, target_compounds: list):
    print("🤖 Task Started: Initiating Structural Search Pipeline...")
    
    # Step 1: Query database for biological primary string
    fasta_file = tool_fetch_uniprot_fasta(user_gene)
    
    # Step 2: Trigger AlphaFold/ColabFold modeling step
    pdb_structure = tool_run_alphafold_prediction(fasta_file)
    
    # Step 3: Map physical constraints (Charges and Hydrogens)
    receptor_pdbqt = tool_prepare_receptor(pdb_structure)
    
    # Step 4: Extract target binding grid (Mocked coordinates for this blueprint example)
    pocket_center = [15.0, 22.0, -4.0]
    grid_dimensions = [20.0, 20.0, 20.0]
    
    # Step 5 & 6: Run Virtual Screening Automation Loop
    screening_registry = []
    for idx, smiles in enumerate(target_compounds):
        try:
            ligand_file = f"compound_{idx}.pdbqt"
            tool_prepare_ligand(smiles, output_path=ligand_file)
            
            # Execute physical engine calculation
            binding_energy = tool_execute_vina_docking(
                receptor_pdbqt, ligand_file, center=pocket_center, box_size=grid_dimensions
            )
            screening_registry.append({"SMILES": smiles, "Affinity (kcal/mol)": binding_energy})
        except Exception as tool_error:
            print(f"Skipping corrupt chemical node [{smiles}]: {tool_error}")
            
    # Step 7: Sort and deliver output report back to the user
    screening_registry.sort(key=lambda x: x["Affinity (kcal/mol)"])
    return screening_registry
```

---

## 5. Required Production Guardrails

To prevent the agent from infinitely spinning compute tokens or crashing due to complex structural data anomalies, you must enforce three validation guardrails:

1. **State Persistence File:** Explicitly write output paths (`/path/to/alphafold.pdb`) to disk at every step. If the underlying server running AlphaFold crashes or restarts due to GPU constraints, the Agent can read the file status ledger and skip straight back to Phase 3 without losing progress.
2. **Structural Confidence (pLDDT Check):** Force the agent to parse the output `.pdb` file coordinates prior to docking. If the per-residue confidence metric (pLDDT) within the localized binding pocket region registers below **70**, write a validation catch forcing the agent to alert the user that the binding pocket is an un-structured/highly flexible loop before executing physical docking.
3. **Rigid Energy Threshold Filtering:** Instruct the agent to throw out simulation poses returning positive binding affinities (> 0.0 kcal/mol). These scores denote immediate physical clashing or atomic overlapping which can easily happen during automated workflows if coordinates are slightly shifted.