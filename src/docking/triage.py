"""
triage.py — Z-score normalisation and consensus hit selection.

Merges Vina and GNINA results, computes Z-scores against a decoy baseline,
flags consensus hits, and saves the ranked hit table.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Z-score computation
# ---------------------------------------------------------------------------

def compute_zscores(scores: pd.Series, mean: float, std: float) -> pd.Series:
    """Compute Z-scores relative to a decoy-set baseline.

    Z = (score − mean) / std

    Parameters
    ----------
    scores : pd.Series
        Series of docking scores (kcal/mol).
    mean : float
        Mean of the decoy score distribution.
    std : float
        Standard deviation of the decoy score distribution.

    Returns
    -------
    pd.Series
        Z-score series with the same index as *scores*.

    Raises
    ------
    ValueError
        If *std* is zero (degenerate distribution).
    """
    if std == 0.0:
        raise ValueError(
            "Standard deviation of decoy scores is zero — cannot compute Z-scores. "
            "Check the decoy set and ensure docking ran correctly."
        )
    return (scores - mean) / std


# ---------------------------------------------------------------------------
# Hit triage
# ---------------------------------------------------------------------------

def triage_hits(
    vina_df: pd.DataFrame,
    gnina_df: pd.DataFrame,
    decoy_mean: float,
    decoy_std: float,
    zscore_cutoff: float = -2.5,
) -> pd.DataFrame:
    """Merge Vina and GNINA results, compute Z-scores, and flag consensus hits.

    Parameters
    ----------
    vina_df : pd.DataFrame
        DataFrame from :func:`~src.docking.screen.run_vina_screen`.
        Required columns: compound_id, smiles, vina_score, pose_file.
    gnina_df : pd.DataFrame
        DataFrame from :func:`~src.docking.gnina.run_gnina_rescore`.
        Required columns: compound_id, vina_score, gnina_score, gnina_cnn_score.
    decoy_mean : float
        Mean decoy Vina score from :func:`~src.docking.screen.run_decoy_screen`.
    decoy_std : float
        Std deviation of decoy Vina scores.
    zscore_cutoff : float
        Z-score threshold for hit flagging (default -2.5). A compound is
        flagged when *both* vina_zscore and gnina_zscore are ≤ this value.
        If all gnina_score values are None, only vina_zscore is used.

    Returns
    -------
    pd.DataFrame
        Columns: compound_id, smiles, vina_score, gnina_score, vina_zscore,
        gnina_zscore, pose_file, consensus_flag.
        Sorted by vina_zscore ascending (most significant hits first).
    """
    # Keep only the columns we need from each frame to avoid merge conflicts
    vina_cols = ["compound_id", "smiles", "vina_score", "pose_file"]
    gnina_cols = ["compound_id", "gnina_score", "gnina_cnn_score"]

    vina_sub = vina_df[vina_cols].copy()

    # gnina_df may contain its own vina_score column — drop it to avoid _x/_y
    gnina_sub = gnina_df[gnina_cols].copy()

    merged = vina_sub.merge(gnina_sub, on="compound_id", how="left")

    # Vina Z-score always uses the decoy baseline
    merged["vina_zscore"] = compute_zscores(
        merged["vina_score"], decoy_mean, decoy_std
    )

    gnina_all_null = merged["gnina_score"].isna().all()

    if gnina_all_null:
        logger.warning(
            "All gnina_score values are None — consensus_flag will use only vina_zscore."
        )
        merged["gnina_zscore"] = float("nan")
        merged["consensus_flag"] = merged["vina_zscore"] <= zscore_cutoff
    else:
        # Use the same decoy baseline for GNINA (both scores are in kcal/mol)
        merged["gnina_zscore"] = compute_zscores(
            merged["gnina_score"].astype(float), decoy_mean, decoy_std
        )
        merged["consensus_flag"] = (
            (merged["vina_zscore"] <= zscore_cutoff)
            & (merged["gnina_zscore"] <= zscore_cutoff)
        )

    # Enforce final column order
    output_columns = [
        "compound_id",
        "smiles",
        "vina_score",
        "gnina_score",
        "vina_zscore",
        "gnina_zscore",
        "pose_file",
        "consensus_flag",
    ]
    merged = merged[output_columns]
    merged.sort_values("vina_zscore", ascending=True, inplace=True)
    merged.reset_index(drop=True, inplace=True)

    n_hits = int(merged["consensus_flag"].sum())
    logger.info(
        "Triage complete: %d/%d compounds are consensus hits (vina_zscore ≤ %.1f "
        "and gnina_zscore ≤ %.1f).",
        n_hits,
        len(merged),
        zscore_cutoff,
        zscore_cutoff,
    )
    return merged


# ---------------------------------------------------------------------------
# Saving results
# ---------------------------------------------------------------------------

def save_hits(df: pd.DataFrame, output_path: str = "results/docking_hits.csv") -> None:
    """Save the triage DataFrame to CSV and print the top 20 hits to stdout.

    Parameters
    ----------
    df : pd.DataFrame
        Full triage DataFrame as returned by :func:`triage_hits`.
    output_path : str
        Destination path for the CSV file (default ``results/docking_hits.csv``).
    """
    out = Path(output_path).resolve()
    os.makedirs(out.parent, exist_ok=True)
    df.to_csv(out, index=False)
    logger.info("Saved %d rows to %s", len(df), out)

    n_consensus = int(df["consensus_flag"].sum()) if "consensus_flag" in df.columns else 0
    print(f"\n{'='*70}")
    print(f"Docking hit triage summary — {len(df)} compounds, {n_consensus} consensus hits")
    print(f"Results saved to: {out}")
    print(f"{'='*70}\n")

    top20 = df.head(20)
    # Pretty-print subset of columns for readability
    display_cols = [c for c in ["compound_id", "vina_score", "gnina_score",
                                 "vina_zscore", "gnina_zscore", "consensus_flag"]
                    if c in top20.columns]
    print("Top 20 compounds by Vina Z-score:")
    print(top20[display_cols].to_string(index=False))
    print()


# ---------------------------------------------------------------------------
# Self-test / usage example
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import tempfile

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    print("=== triage.py self-test ===")

    # Build synthetic screening results
    vina_df = pd.DataFrame(
        {
            "compound_id": [f"cpd_{i:03d}" for i in range(1, 11)],
            "smiles": ["CCO"] * 10,
            "vina_score": [-9.5, -8.8, -8.1, -7.5, -7.0, -6.5, -6.0, -5.5, -5.0, -4.5],
            "pose_file": [f"/tmp/cpd_{i:03d}_pose.pdbqt" for i in range(1, 11)],
        }
    )

    gnina_df = pd.DataFrame(
        {
            "compound_id": [f"cpd_{i:03d}" for i in range(1, 11)],
            "vina_score": vina_df["vina_score"].tolist(),
            "gnina_score": [-10.1, -9.2, -8.5, -7.8, -7.1, -6.4, -5.9, -5.2, -4.8, -4.2],
            "gnina_cnn_score": [0.85, 0.72, 0.68, 0.61, 0.55, 0.48, 0.43, 0.38, 0.31, 0.27],
        }
    )

    # Decoy baseline: mean=-6.0, std=1.0
    decoy_mean = -6.0
    decoy_std = 1.0

    result = triage_hits(vina_df, gnina_df, decoy_mean, decoy_std, zscore_cutoff=-2.5)

    print("triage_hits output columns:", list(result.columns))
    print(result.to_string(index=False))
    print(f"\nConsensus hits: {result['consensus_flag'].sum()}")

    # Test save_hits
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = os.path.join(tmpdir, "docking_hits.csv")
        save_hits(result, output_path=csv_path)
        print(f"CSV saved, {os.path.getsize(csv_path)} bytes")

    # Test gnina_all_null fallback
    gnina_null_df = gnina_df.copy()
    gnina_null_df["gnina_score"] = None
    gnina_null_df["gnina_cnn_score"] = None
    result_no_gnina = triage_hits(
        vina_df, gnina_null_df, decoy_mean, decoy_std, zscore_cutoff=-2.5
    )
    print("\ngno-GNINA fallback consensus_flag (vina-only):")
    print(result_no_gnina[["compound_id", "vina_zscore", "consensus_flag"]].to_string(index=False))

    print("\n=== self-test complete ===")
