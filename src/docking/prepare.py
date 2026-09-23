"""
prepare.py — Receptor and ligand preparation for AutoDock Vina screening.

Converts AF3 CIF output → PDBQT receptor, and SMILES → PDBQT ligands.
"""

from __future__ import annotations

import logging
import os
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Receptor preparation
# ---------------------------------------------------------------------------

def extract_chain(cif_path: str, chain_id: str, out_pdb_path: str) -> str:
    """Extract a single chain from an AF3 CIF file and write it as PDB.

    Parameters
    ----------
    cif_path : str
        Path to the input CIF file produced by AlphaFold 3.
    chain_id : str
        Chain identifier to extract (e.g. "A").
    out_pdb_path : str
        Destination path for the output PDB file.

    Returns
    -------
    str
        Absolute path to the written PDB file.
    """
    from Bio.PDB import MMCIFParser, PDBIO, Select

    class ChainSelect(Select):
        def __init__(self, chain: str) -> None:
            self.chain = chain

        def accept_chain(self, chain) -> bool:  # type: ignore[override]
            return chain.get_id() == self.chain

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("receptor", cif_path)

    out_pdb_path = str(Path(out_pdb_path).resolve())
    os.makedirs(os.path.dirname(out_pdb_path) or ".", exist_ok=True)

    io = PDBIO()
    io.set_structure(structure)
    io.save(out_pdb_path, ChainSelect(chain_id))

    logger.info("Extracted chain %s from %s → %s", chain_id, cif_path, out_pdb_path)
    return out_pdb_path


def fix_receptor(pdb_path: str, out_pdb_path: str) -> str:
    """Add missing atoms / hydrogens and remove heteroatoms via PDBFixer.

    Adds hydrogens at pH 7.4, finds missing residues / heavy atoms,
    and removes all HETATM records (non-protein residues).

    Parameters
    ----------
    pdb_path : str
        Path to the input PDB file.
    out_pdb_path : str
        Destination path for the fixed PDB file.

    Returns
    -------
    str
        Absolute path to the fixed PDB file.
    """
    from pdbfixer import PDBFixer
    import openmm.app as app

    fixer = PDBFixer(filename=pdb_path)

    # Find and fill missing residues / atoms
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)

    out_pdb_path = str(Path(out_pdb_path).resolve())
    os.makedirs(os.path.dirname(out_pdb_path) or ".", exist_ok=True)

    with open(out_pdb_path, "w") as fh:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, fh)

    logger.info("Fixed receptor %s → %s", pdb_path, out_pdb_path)
    return out_pdb_path


def receptor_to_pdbqt(pdb_path: str, out_pdbqt_path: str) -> str:
    """Convert a fixed PDB receptor to PDBQT format using Meeko.

    Parameters
    ----------
    pdb_path : str
        Path to the cleaned/fixed PDB file.
    out_pdbqt_path : str
        Destination path for the PDBQT file.

    Returns
    -------
    str
        Absolute path to the written PDBQT file.
    """
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from rdkit import Chem

    mol = Chem.MolFromPDBFile(pdb_path, removeHs=False, sanitize=True)
    if mol is None:
        raise ValueError(f"RDKit could not parse PDB file: {pdb_path}")

    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)

    out_pdbqt_path = str(Path(out_pdbqt_path).resolve())
    os.makedirs(os.path.dirname(out_pdbqt_path) or ".", exist_ok=True)

    pdbqt_string, _, _ = PDBQTWriterLegacy.write_string(mol_setups[0])
    with open(out_pdbqt_path, "w") as fh:
        fh.write(pdbqt_string)

    logger.info("Converted receptor %s → %s", pdb_path, out_pdbqt_path)
    return out_pdbqt_path


def prepare_receptor(cif_path: str, chain_id: str, output_dir: str) -> str:
    """Orchestrate CIF → PDB → fixed PDB → PDBQT pipeline for a receptor.

    Parameters
    ----------
    cif_path : str
        Path to the AF3 CIF output file.
    chain_id : str
        Chain ID of the protein receptor (e.g. "A").
    output_dir : str
        Directory where intermediate and final files are written.

    Returns
    -------
    str
        Absolute path to the final receptor PDBQT file.
    """
    output_dir = str(Path(output_dir).resolve())
    os.makedirs(output_dir, exist_ok=True)

    stem = Path(cif_path).stem
    raw_pdb = os.path.join(output_dir, f"{stem}_chain{chain_id}.pdb")
    fixed_pdb = os.path.join(output_dir, f"{stem}_chain{chain_id}_fixed.pdb")
    receptor_pdbqt = os.path.join(output_dir, f"{stem}_chain{chain_id}.pdbqt")

    extract_chain(cif_path, chain_id, raw_pdb)
    fix_receptor(raw_pdb, fixed_pdb)
    receptor_to_pdbqt(fixed_pdb, receptor_pdbqt)

    return receptor_pdbqt


# ---------------------------------------------------------------------------
# Ligand preparation
# ---------------------------------------------------------------------------

def smiles_to_pdbqt(smiles: str, compound_id: str, output_dir: str) -> Optional[str]:
    """Convert a SMILES string to an energy-minimized PDBQT ligand file.

    Uses RDKit ETKDGv3 embedding + MMFF94 minimisation followed by Meeko
    to write the PDBQT.  Returns None on any failure (logs a warning).

    Parameters
    ----------
    smiles : str
        SMILES string for the compound.
    compound_id : str
        Unique identifier used to name the output file.
    output_dir : str
        Directory where the PDBQT file is written.

    Returns
    -------
    str or None
        Absolute path to the PDBQT file, or None if preparation failed.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from meeko import MoleculePreparation, PDBQTWriterLegacy

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning("Could not parse SMILES for %s: %s", compound_id, smiles)
            return None

        mol = Chem.AddHs(mol)

        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        result = AllChem.EmbedMolecule(mol, params)
        if result == -1:
            logger.warning("3-D embedding failed for %s (SMILES: %s)", compound_id, smiles)
            return None

        ff_result = AllChem.MMFFOptimizeMolecule(mol, forceField="MMFF94")
        if ff_result not in (0, 1):
            logger.warning("MMFF optimisation did not converge for %s", compound_id)
            # Non-fatal — proceed with the embedded geometry

        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)

        output_dir = str(Path(output_dir).resolve())
        os.makedirs(output_dir, exist_ok=True)
        out_path = os.path.join(output_dir, f"{compound_id}.pdbqt")

        pdbqt_string, _, _ = PDBQTWriterLegacy.write_string(mol_setups[0])
        with open(out_path, "w") as fh:
            fh.write(pdbqt_string)

        return out_path

    except Exception as exc:  # noqa: BLE001
        logger.warning("smiles_to_pdbqt failed for %s: %s", compound_id, exc)
        return None


# Worker function at module level so ProcessPoolExecutor can pickle it
def _prepare_one(args: tuple[str, str, str]) -> Optional[str]:
    smiles, compound_id, output_dir = args
    return smiles_to_pdbqt(smiles, compound_id, output_dir)


def prepare_library(
    smi_file: str,
    output_dir: str,
    n_workers: int = 4,
) -> list[str]:
    """Prepare a compound library from a SMILES file in parallel.

    Each non-comment line should be formatted as::

        SMILES  compound_id

    Blank lines and lines starting with ``#`` are skipped.

    Parameters
    ----------
    smi_file : str
        Path to the SMILES file.
    output_dir : str
        Directory where PDBQT files are written.
    n_workers : int
        Number of worker processes (default 4).

    Returns
    -------
    list[str]
        Sorted list of absolute paths to successfully prepared PDBQT files.
    """
    entries: list[tuple[str, str, str]] = []
    with open(smi_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                logger.warning("Skipping malformed SMILES line: %r", line)
                continue
            smiles, compound_id = parts[0], parts[1]
            entries.append((smiles, compound_id, output_dir))

    logger.info("Preparing %d compounds with %d workers…", len(entries), n_workers)
    successful: list[str] = []

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        future_to_id = {
            executor.submit(_prepare_one, entry): entry[1] for entry in entries
        }
        for future in as_completed(future_to_id):
            compound_id = future_to_id[future]
            try:
                result = future.result()
            except Exception as exc:  # noqa: BLE001
                logger.warning("Worker failed for %s: %s", compound_id, exc)
                result = None

            if result is not None:
                successful.append(result)

    logger.info(
        "Library preparation complete: %d/%d compounds succeeded.",
        len(successful),
        len(entries),
    )
    return sorted(successful)


# ---------------------------------------------------------------------------
# Self-test / usage example
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import tempfile

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    print("=== prepare.py self-test ===")

    with tempfile.TemporaryDirectory() as tmpdir:
        # Test smiles_to_pdbqt with a simple drug-like molecule (ibuprofen)
        test_smiles = "CC(Cc1ccc(cc1)C(C)C(=O)O)C"
        test_id = "ibuprofen_test"
        result = smiles_to_pdbqt(test_smiles, test_id, tmpdir)
        if result and os.path.exists(result):
            print(f"smiles_to_pdbqt OK → {result}")
        else:
            print("smiles_to_pdbqt FAILED (meeko/rdkit may not be installed in this env)")

        # Test prepare_library with an inline SMILES file
        smi_content = (
            "CC(=O)Oc1ccccc1C(=O)O aspirin\n"
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C caffeine\n"
        )
        smi_path = os.path.join(tmpdir, "test.smi")
        with open(smi_path, "w") as fh:
            fh.write(smi_content)

        paths = prepare_library(smi_path, os.path.join(tmpdir, "ligands"), n_workers=2)
        print(f"prepare_library returned {len(paths)} path(s): {paths}")

    print("=== self-test complete ===")
