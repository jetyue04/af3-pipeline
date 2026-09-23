"""
gnina.py — GNINA rescoring for top AutoDock Vina hits.

Uses the GNINA binary (must be on PATH) to apply CNN-based rescoring on
top-ranked poses from the Vina screening pass.
"""

from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

# Regex patterns for GNINA stdout parsing
_AFFINITY_RE = re.compile(r"Affinity:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)")
_CNN_SCORE_RE = re.compile(r"CNNscore:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)")


# ---------------------------------------------------------------------------
# Availability check
# ---------------------------------------------------------------------------

def check_gnina() -> bool:
    """Return True if the ``gnina`` binary is available on PATH.

    Logs a warning when GNINA is not found so callers can fall back
    gracefully.

    Returns
    -------
    bool
        ``True`` if gnina is on PATH, ``False`` otherwise.
    """
    found = shutil.which("gnina") is not None
    if not found:
        logger.warning(
            "gnina binary not found on PATH.  CNN rescoring will be skipped. "
            "Install GNINA from https://github.com/gnina/gnina/releases and "
            "ensure it is executable and on your PATH."
        )
    return found


# ---------------------------------------------------------------------------
# Single-pose rescoring
# ---------------------------------------------------------------------------

def rescore_with_gnina(
    receptor_pdbqt: str,
    pose_file: str,
    out_dir: str,
) -> dict:
    """Rescore a single docked pose with GNINA (score_only mode).

    Runs::

        gnina -r receptor.pdbqt -l pose.pdbqt --score_only

    and parses the Affinity and CNNscore values from stdout.

    Parameters
    ----------
    receptor_pdbqt : str
        Path to the receptor PDBQT file.
    pose_file : str
        Path to the pose PDBQT file (typically the best Vina pose).
    out_dir : str
        Directory for any incidental output files.  Not used in
        ``--score_only`` mode but retained for API consistency.

    Returns
    -------
    dict
        ``{"gnina_score": float, "gnina_cnn_score": float}``

    Raises
    ------
    RuntimeError
        If the GNINA process returns a non-zero exit code or the expected
        score values cannot be parsed from its output.
    """
    os.makedirs(out_dir, exist_ok=True)

    cmd = [
        "gnina",
        "-r", receptor_pdbqt,
        "-l", pose_file,
        "--score_only",
    ]

    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"GNINA timed out rescoring {pose_file}"
        ) from exc
    except FileNotFoundError as exc:
        raise RuntimeError(
            "gnina binary not found — install GNINA and add it to PATH."
        ) from exc

    if proc.returncode != 0:
        raise RuntimeError(
            f"GNINA exited with code {proc.returncode}.\n"
            f"stderr: {proc.stderr.strip()}"
        )

    stdout = proc.stdout

    affinity_match = _AFFINITY_RE.search(stdout)
    cnn_match = _CNN_SCORE_RE.search(stdout)

    if affinity_match is None or cnn_match is None:
        raise RuntimeError(
            f"Could not parse GNINA scores from output:\n{stdout[:500]}"
        )

    return {
        "gnina_score": float(affinity_match.group(1)),
        "gnina_cnn_score": float(cnn_match.group(1)),
    }


# ---------------------------------------------------------------------------
# Batch rescoring
# ---------------------------------------------------------------------------

def run_gnina_rescore(
    receptor_pdbqt: str,
    vina_df: pd.DataFrame,
    out_dir: str,
    top_n: int = 500,
) -> pd.DataFrame:
    """Rescore the top-N Vina hits with GNINA.

    Takes the best ``top_n`` compounds by Vina score, rescores each with
    GNINA, and returns a merged DataFrame.  If GNINA is unavailable, the
    function returns an augmented copy of ``vina_df`` (top-N rows) with
    ``gnina_score`` and ``gnina_cnn_score`` set to ``None`` and logs a
    warning.

    Parameters
    ----------
    receptor_pdbqt : str
        Path to the receptor PDBQT file.
    vina_df : pd.DataFrame
        DataFrame as returned by :func:`~src.docking.screen.run_vina_screen`.
        Must contain columns: compound_id, vina_score, pose_file.
    out_dir : str
        Directory where GNINA output files are written.
    top_n : int
        Number of top Vina hits to rescore (default 500).

    Returns
    -------
    pd.DataFrame
        Columns: compound_id, vina_score, gnina_score, gnina_cnn_score.
        Sorted by vina_score ascending.
    """
    os.makedirs(out_dir, exist_ok=True)

    # Sort by Vina score (most favourable = most negative first) and take top N
    top_df = (
        vina_df.sort_values("vina_score", ascending=True)
        .head(top_n)
        .reset_index(drop=True)
    )

    if not check_gnina():
        logger.warning(
            "GNINA unavailable — returning top-%d Vina results with gnina_score=None.",
            len(top_df),
        )
        result = top_df[["compound_id", "vina_score"]].copy()
        result["gnina_score"] = None
        result["gnina_cnn_score"] = None
        return result[["compound_id", "vina_score", "gnina_score", "gnina_cnn_score"]]

    rows: list[dict] = []
    for _, row in top_df.iterrows():
        compound_id = row["compound_id"]
        vina_score = row["vina_score"]
        pose_file = row["pose_file"]

        try:
            gnina_result = rescore_with_gnina(receptor_pdbqt, pose_file, out_dir)
            rows.append(
                {
                    "compound_id": compound_id,
                    "vina_score": vina_score,
                    "gnina_score": gnina_result["gnina_score"],
                    "gnina_cnn_score": gnina_result["gnina_cnn_score"],
                }
            )
        except Exception as exc:  # noqa: BLE001
            logger.warning("GNINA rescoring failed for %s: %s", compound_id, exc)
            rows.append(
                {
                    "compound_id": compound_id,
                    "vina_score": vina_score,
                    "gnina_score": None,
                    "gnina_cnn_score": None,
                }
            )

    result_df = pd.DataFrame(
        rows, columns=["compound_id", "vina_score", "gnina_score", "gnina_cnn_score"]
    )
    result_df.sort_values("vina_score", ascending=True, inplace=True)
    result_df.reset_index(drop=True, inplace=True)

    logger.info(
        "GNINA rescore complete: %d/%d compounds have valid GNINA scores.",
        result_df["gnina_score"].notna().sum(),
        len(result_df),
    )
    return result_df


# ---------------------------------------------------------------------------
# Self-test / usage example
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    print("=== gnina.py self-test ===")

    gnina_available = check_gnina()
    print(f"GNINA on PATH: {gnina_available}")

    # Demonstrate output contract
    dummy_vina_df = pd.DataFrame(
        {
            "compound_id": ["cpd_001", "cpd_002", "cpd_003"],
            "smiles": ["CC(=O)O", "c1ccccc1", "CCO"],
            "vina_score": [-8.5, -7.2, -6.1],
            "pose_file": ["/tmp/cpd_001_pose.pdbqt", "/tmp/cpd_002_pose.pdbqt", "/tmp/cpd_003_pose.pdbqt"],
        }
    )

    print("\nDummy vina_df:")
    print(dummy_vina_df.to_string(index=False))

    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        # This will trigger the "gnina unavailable" path if gnina is not installed
        result_df = run_gnina_rescore(
            receptor_pdbqt="/tmp/receptor.pdbqt",
            vina_df=dummy_vina_df,
            out_dir=tmpdir,
            top_n=3,
        )

    print("\nrun_gnina_rescore output columns:", list(result_df.columns))
    print(result_df.to_string(index=False))
    print("=== self-test complete ===")
