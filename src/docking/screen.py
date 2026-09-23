"""
screen.py — AutoDock Vina batch virtual screening.

Defines the docking grid, runs individual and batch Vina jobs, and provides
a decoy-screen helper for Z-score baseline calibration.
"""

from __future__ import annotations

import logging
import os
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Grid definition
# ---------------------------------------------------------------------------

def define_grid(
    contact_residues: list[int],
    cif_path: str,
    chain_id: str,
    padding: float = 4.0,
) -> dict:
    """Compute a Vina search box centred on interface contact residues.

    The centroid of Cα atoms for the given residue numbers is computed and
    a cubic box large enough to enclose those atoms (plus padding on each
    side) is returned.

    Parameters
    ----------
    contact_residues : list[int]
        Sequence numbers of the interface contact residues.
    cif_path : str
        Path to the AF3 CIF file containing the receptor structure.
    chain_id : str
        Chain containing the contact residues.
    padding : float
        Extra space (Å) added on every side of the bounding box.

    Returns
    -------
    dict
        Keys: center_x, center_y, center_z, size_x, size_y, size_z (floats).
    """
    import numpy as np
    from Bio.PDB import MMCIFParser

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("grid_calc", cif_path)

    ca_coords: list[tuple[float, float, float]] = []
    for model in structure:
        try:
            chain = model[chain_id]
        except KeyError:
            continue
        for residue in chain:
            if residue.get_id()[1] in contact_residues:
                if "CA" in residue:
                    coord = residue["CA"].get_vector().get_array()
                    ca_coords.append(tuple(coord))

    if not ca_coords:
        raise ValueError(
            f"No Cα atoms found for residues {contact_residues} in chain "
            f"{chain_id} of {cif_path}."
        )

    coords = np.array(ca_coords)  # shape (N, 3)
    centroid = coords.mean(axis=0)
    span = coords.max(axis=0) - coords.min(axis=0)
    box_size = span + 2 * padding  # symmetric padding on both sides

    return {
        "center_x": float(centroid[0]),
        "center_y": float(centroid[1]),
        "center_z": float(centroid[2]),
        "size_x": float(box_size[0]),
        "size_y": float(box_size[1]),
        "size_z": float(box_size[2]),
    }


# ---------------------------------------------------------------------------
# Single-ligand Vina docking
# ---------------------------------------------------------------------------

def run_vina(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    grid: dict,
    out_dir: str,
    exhaustiveness: int = 8,
) -> dict:
    """Dock one ligand against the receptor and return the best pose.

    Parameters
    ----------
    receptor_pdbqt : str
        Path to the receptor PDBQT file.
    ligand_pdbqt : str
        Path to the ligand PDBQT file.
    grid : dict
        Grid dictionary as returned by :func:`define_grid`.
    out_dir : str
        Directory where the pose file is written.
    exhaustiveness : int
        Vina exhaustiveness parameter (8 for screening, 32 for rescoring).

    Returns
    -------
    dict
        ``{"score": float, "pose_file": str}`` — Vina binding energy (kcal/mol)
        and absolute path to the written pose PDBQT.

    Raises
    ------
    RuntimeError
        Propagated if Vina itself raises an unrecoverable error.
    """
    from vina import Vina

    os.makedirs(out_dir, exist_ok=True)
    compound_id = Path(ligand_pdbqt).stem
    pose_file = os.path.join(str(Path(out_dir).resolve()), f"{compound_id}_pose.pdbqt")

    v = Vina(sf_name="vina", cpu=1, seed=42, verbosity=0)
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(
        center=[grid["center_x"], grid["center_y"], grid["center_z"]],
        box_size=[grid["size_x"], grid["size_y"], grid["size_z"]],
    )
    v.dock(exhaustiveness=exhaustiveness, n_poses=5)
    best_score = float(v.energies(n_poses=1)[0][0])
    v.write_poses(pose_file, n_poses=1, overwrite=True)

    return {"score": best_score, "pose_file": pose_file}


# ---------------------------------------------------------------------------
# Parallel batch screening
# ---------------------------------------------------------------------------

# Module-level worker — must be picklable for ProcessPoolExecutor
def _vina_worker(
    args: tuple[str, str, str, dict, str, int],
) -> Optional[dict]:
    receptor_pdbqt, ligand_pdbqt, compound_id, grid, out_dir, exhaustiveness = args
    try:
        result = run_vina(receptor_pdbqt, ligand_pdbqt, grid, out_dir, exhaustiveness)
        return {
            "compound_id": compound_id,
            "score": result["score"],
            "pose_file": result["pose_file"],
        }
    except Exception as exc:  # noqa: BLE001
        logger.warning("Vina failed for %s: %s", compound_id, exc)
        return None


def run_vina_screen(
    receptor_pdbqt: str,
    ligand_pdbqt_list: list[str],
    grid: dict,
    out_dir: str,
    n_workers: int = 4,
) -> pd.DataFrame:
    """Run AutoDock Vina on a list of ligands in parallel.

    Failures for individual compounds are logged and skipped; the screen
    never crashes due to a single bad ligand.

    The ``smiles`` column is populated as an empty string here because
    SMILES are not stored inside PDBQT files.  Callers that need SMILES
    should merge the returned DataFrame with their compound registry.

    Parameters
    ----------
    receptor_pdbqt : str
        Path to the receptor PDBQT file.
    ligand_pdbqt_list : list[str]
        Paths to ligand PDBQT files.
    grid : dict
        Grid dictionary as returned by :func:`define_grid`.
    out_dir : str
        Directory where pose files are written.
    n_workers : int
        Number of parallel worker processes.

    Returns
    -------
    pd.DataFrame
        Columns: compound_id, smiles, vina_score, pose_file.
        Sorted by vina_score ascending (most favourable first).
    """
    os.makedirs(out_dir, exist_ok=True)

    work_items = [
        (
            receptor_pdbqt,
            lig,
            Path(lig).stem,
            grid,
            out_dir,
            8,  # exhaustiveness for screening pass
        )
        for lig in ligand_pdbqt_list
    ]

    logger.info(
        "Starting Vina screen: %d ligands, %d workers", len(work_items), n_workers
    )

    rows: list[dict] = []
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        future_to_id = {
            executor.submit(_vina_worker, item): item[2] for item in work_items
        }
        for future in as_completed(future_to_id):
            compound_id = future_to_id[future]
            try:
                result = future.result()
            except Exception as exc:  # noqa: BLE001
                logger.warning("Worker raised for %s: %s", compound_id, exc)
                result = None

            if result is not None:
                rows.append(
                    {
                        "compound_id": result["compound_id"],
                        "smiles": "",  # populated by caller from compound registry
                        "vina_score": result["score"],
                        "pose_file": result["pose_file"],
                    }
                )

    df = pd.DataFrame(rows, columns=["compound_id", "smiles", "vina_score", "pose_file"])
    df.sort_values("vina_score", ascending=True, inplace=True)
    df.reset_index(drop=True, inplace=True)

    logger.info(
        "Screen complete: %d/%d ligands scored.", len(df), len(ligand_pdbqt_list)
    )
    return df


# ---------------------------------------------------------------------------
# Decoy screen for Z-score calibration
# ---------------------------------------------------------------------------

def run_decoy_screen(
    receptor_pdbqt: str,
    decoy_smi_file: str,
    grid: dict,
    out_dir: str,
) -> tuple[float, float]:
    """Screen the decoy set and return mean/std binding scores for Z-score calibration.

    Prepares each decoy SMILES on-the-fly and runs single-threaded Vina
    docking (the decoy set is small enough that parallelism adds overhead).

    Parameters
    ----------
    receptor_pdbqt : str
        Path to the receptor PDBQT file.
    decoy_smi_file : str
        Path to the decoy SMILES file (same "SMILES id" format).
    grid : dict
        Grid dictionary as returned by :func:`define_grid`.
    out_dir : str
        Directory for temporary decoy PDBQT and pose files.

    Returns
    -------
    tuple[float, float]
        ``(mean_score, std_score)`` of all successful decoy docking scores.

    Raises
    ------
    RuntimeError
        If no decoys could be scored.
    """
    import numpy as np
    from src.docking.prepare import smiles_to_pdbqt

    os.makedirs(out_dir, exist_ok=True)
    decoy_pdbqt_dir = os.path.join(out_dir, "decoy_pdbqt")
    decoy_pose_dir = os.path.join(out_dir, "decoy_poses")
    os.makedirs(decoy_pdbqt_dir, exist_ok=True)
    os.makedirs(decoy_pose_dir, exist_ok=True)

    decoys: list[tuple[str, str]] = []
    with open(decoy_smi_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            decoys.append((parts[0], parts[1]))

    scores: list[float] = []
    for smiles, compound_id in decoys:
        try:
            lig_path = smiles_to_pdbqt(smiles, compound_id, decoy_pdbqt_dir)
            if lig_path is None:
                continue
            result = run_vina(receptor_pdbqt, lig_path, grid, decoy_pose_dir, exhaustiveness=8)
            scores.append(result["score"])
        except Exception as exc:  # noqa: BLE001
            logger.warning("Decoy docking failed for %s: %s", compound_id, exc)

    if not scores:
        raise RuntimeError("No decoys could be scored — cannot calibrate Z-scores.")

    mean_score = float(np.mean(scores))
    std_score = float(np.std(scores))
    logger.info(
        "Decoy screen: n=%d, mean=%.3f, std=%.3f", len(scores), mean_score, std_score
    )
    return mean_score, std_score


# ---------------------------------------------------------------------------
# Self-test / usage example
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import json

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    print("=== screen.py self-test ===")

    # Demonstrate define_grid signature (requires a real CIF file to run fully)
    print("define_grid signature: define_grid(contact_residues, cif_path, chain_id, padding=4.0)")
    print("Returns: {center_x, center_y, center_z, size_x, size_y, size_z}")

    # Demonstrate run_vina_screen signature
    print("\nrun_vina_screen returns a DataFrame with columns:")
    empty_df = pd.DataFrame(columns=["compound_id", "smiles", "vina_score", "pose_file"])
    print(empty_df.dtypes.to_string())

    print("\nrun_decoy_screen returns (mean_score, std_score) tuple for Z-score calibration.")
    print("=== self-test complete ===")
