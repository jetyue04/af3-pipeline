"""Wrapper around prodigy-prot for binding ΔG / Kd prediction from AF3 CIF files."""
import copy
import logging
import os
import warnings
from pathlib import Path
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)

# Try to import prodigy_prot once at module load so callers can check availability
try:
    import prodigy_prot  # type: ignore  # noqa: F401
    _PRODIGY_AVAILABLE = True
except ImportError:
    _PRODIGY_AVAILABLE = False
    logger.warning(
        "prodigy-prot is not installed. run_prodigy() will return None values. "
        "Install with: pip install prodigy-prot"
    )


def _extract_protein_model(cif_path: str, chain_ids: set[str]):
    """Parse a CIF file and return a Model containing only the target chains.

    Non-protein chains (ligands, ions, water) are dropped so they do not
    pollute PRODIGY's SASA/contact calculation.
    """
    from Bio.PDB import MMCIFParser, Model, Structure  # type: ignore

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("struct", cif_path)
    src_model = next(iter(structure))

    new_structure = Structure.Structure("prodigy_sel")
    model = Model.Model(0)
    new_structure.add(model)
    for chain in src_model.get_chains():
        if chain.get_id() in chain_ids:
            model.add(copy.deepcopy(chain))
    return model


def run_prodigy(
    cif_path: str,
    chain_a_id: str,
    chain_b_id: str,
) -> dict:
    """Predict binding ΔG and Kd for an interface using prodigy-prot.

    Parameters
    ----------
    cif_path:
        Path to the AF3 CIF model file.
    chain_a_id:
        Chain identifier for chain A.
    chain_b_id:
        Chain identifier for chain B.

    Returns
    -------
    dict with keys:
        dg_kcal_mol (float | None): binding free energy in kcal/mol
        kd_molar   (float | None): dissociation constant in molar units
    """
    if not _PRODIGY_AVAILABLE:
        warnings.warn(
            "prodigy-prot is not installed; returning None values.",
            RuntimeWarning,
            stacklevel=2,
        )
        return {"dg_kcal_mol": None, "kd_molar": None}

    try:
        model = _extract_protein_model(cif_path, {chain_a_id, chain_b_id})

        from prodigy_prot.modules.prodigy import Prodigy  # type: ignore

        prodigy = Prodigy(model, selection=[chain_a_id, chain_b_id])
        prodigy.predict(temp=25.0)

        return {
            "dg_kcal_mol": float(prodigy.ba_val),
            "kd_molar": float(prodigy.kd_val),
        }

    except Exception as exc:
        logger.error(
            "prodigy prediction failed for %s (%s/%s): %s",
            cif_path, chain_a_id, chain_b_id, exc,
        )
        return {"dg_kcal_mol": None, "kd_molar": None}


def run_prodigy_batch(jobs: list[dict]) -> pd.DataFrame:
    """Run prodigy on a list of interface jobs.

    Parameters
    ----------
    jobs:
        List of dicts, each with keys:
            job_name (str): identifier for the job
            cif_path (str): path to the AF3 CIF model file
            chain_a  (str): chain A identifier
            chain_b  (str): chain B identifier

    Returns
    -------
    pd.DataFrame with columns: job_name, chain_a, chain_b, dg_kcal_mol, kd_molar
    """
    rows: list[dict] = []
    for job in jobs:
        job_name = job["job_name"]
        cif_path = job["cif_path"]
        chain_a = job["chain_a"]
        chain_b = job["chain_b"]

        logger.info("Running prodigy: %s (%s / %s)", job_name, chain_a, chain_b)
        result = run_prodigy(cif_path, chain_a, chain_b)
        rows.append(
            {
                "job_name": job_name,
                "chain_a": chain_a,
                "chain_b": chain_b,
                "dg_kcal_mol": result["dg_kcal_mol"],
                "kd_molar": result["kd_molar"],
            }
        )

    return pd.DataFrame(
        rows,
        columns=["job_name", "chain_a", "chain_b", "dg_kcal_mol", "kd_molar"],
    )


if __name__ == "__main__":
    import sys

    logging.basicConfig(level=logging.INFO, stream=sys.stderr)

    output_dir = os.environ.get(
        "AF3_OUTPUT_DIR",
        str(Path(__file__).resolve().parent.parent.parent / "output"),
    )

    cif_files = list(Path(output_dir).glob("**/*_model.cif"))
    if not cif_files:
        # Broader search for any cif
        cif_files = list(Path(output_dir).glob("**/*.cif"))

    if not cif_files:
        print("No CIF files found for self-test — skipping.")
        sys.exit(0)

    test_cif = str(cif_files[0])
    print(f"Testing run_prodigy with: {test_cif}")

    # Detect chain IDs from file
    from Bio.PDB import MMCIFParser  # type: ignore
    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure("s", test_cif)
    model = next(iter(struct))
    chain_ids = [c.get_id() for c in model.get_chains()]
    print(f"Chains: {chain_ids}")

    if len(chain_ids) >= 2:
        result = run_prodigy(test_cif, chain_ids[0], chain_ids[1])
        print("run_prodigy result:", result)

        jobs = [
            {
                "job_name": Path(test_cif).stem,
                "cif_path": test_cif,
                "chain_a": chain_ids[0],
                "chain_b": chain_ids[1],
            }
        ]
        df = run_prodigy_batch(jobs)
        print("\nrun_prodigy_batch result:")
        print(df.to_string())
    else:
        print("Only one chain found; skipping prodigy test.")
