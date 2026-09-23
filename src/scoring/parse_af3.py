"""Parse AlphaFold 3 summary_confidences.json files into a score table."""
import os
import json
import glob
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def infer_nucleotide(job_name: str) -> str:
    """Infer nucleotide state from job name: 'GTP', 'GDP', or 'none'."""
    lower = job_name.lower()
    if "gtp" in lower:
        return "GTP"
    if "gdp" in lower:
        return "GDP"
    return "none"


def load_job(job_dir: str) -> dict:
    """Load summary_confidences.json from a single job directory.

    Parameters
    ----------
    job_dir:
        Path to an AF3 job output directory (e.g. output/my_job/).

    Returns
    -------
    dict with keys:
        job_name, iptm, ptm, ranking_score, fraction_disordered, has_clash,
        chain_iptm_mean, cross_chain_pae_min, nucleotide
    """
    job_path = Path(job_dir)
    job_name = job_path.name

    # AF3 places summary_confidences.json directly inside the job directory
    summary_path = job_path / f"{job_name}_summary_confidences.json"
    if not summary_path.exists():
        # Fallback: search for any *_summary_confidences.json
        candidates = list(job_path.glob("*_summary_confidences.json"))
        if not candidates:
            raise FileNotFoundError(
                f"No summary_confidences.json found in {job_dir}"
            )
        summary_path = candidates[0]

    with open(summary_path, "r") as fh:
        data = json.load(fh)

    # chain_iptm_mean: mean of non-null values in chain_iptm list
    chain_iptm_raw = data.get("chain_iptm", [])
    non_null_iptm = [v for v in chain_iptm_raw if v is not None]
    chain_iptm_mean = float(np.mean(non_null_iptm)) if non_null_iptm else float("nan")

    # cross_chain_pae_min: minimum off-diagonal value in chain_pair_pae_min
    pae_matrix = data.get("chain_pair_pae_min", [])
    cross_chain_pae_min = float("nan")
    if pae_matrix:
        n = len(pae_matrix)
        off_diag = [
            pae_matrix[i][j]
            for i in range(n)
            for j in range(n)
            if i != j and pae_matrix[i][j] is not None
        ]
        if off_diag:
            cross_chain_pae_min = float(min(off_diag))

    return {
        "job_name": job_name,
        "iptm": data.get("iptm"),
        "ptm": data.get("ptm"),
        "ranking_score": data.get("ranking_score"),
        "fraction_disordered": data.get("fraction_disordered"),
        "has_clash": data.get("has_clash"),
        "chain_iptm_mean": chain_iptm_mean,
        "cross_chain_pae_min": cross_chain_pae_min,
        "nucleotide": infer_nucleotide(job_name),
    }


def parse_af3(output_dir: str | None = None) -> pd.DataFrame:
    """Scan output_dir for AF3 job subdirectories and return a score DataFrame.

    Precedence for output directory:
    1. ``output_dir`` argument (if not None)
    2. ``AF3_OUTPUT_DIR`` environment variable
    3. ``../output`` relative to this file's location

    Parameters
    ----------
    output_dir:
        Optional explicit path to the AF3 output root.

    Returns
    -------
    pd.DataFrame with columns:
        job_name, iptm, ptm, ranking_score, fraction_disordered, has_clash,
        chain_iptm_mean, cross_chain_pae_min, nucleotide
    """
    if output_dir is None:
        output_dir = os.environ.get(
            "AF3_OUTPUT_DIR",
            str(Path(__file__).resolve().parent.parent.parent / "output"),
        )

    output_path = Path(output_dir)
    if not output_path.exists():
        logger.warning("AF3 output directory does not exist: %s", output_path)
        return _empty_dataframe()

    # Each subdirectory is one job
    job_dirs = [p for p in sorted(output_path.iterdir()) if p.is_dir()]
    if not job_dirs:
        logger.warning("No job subdirectories found in %s", output_path)
        return _empty_dataframe()

    rows: list[dict] = []
    for job_dir in job_dirs:
        try:
            row = load_job(str(job_dir))
            rows.append(row)
        except FileNotFoundError as exc:
            logger.warning("Skipping %s: %s", job_dir.name, exc)
        except (KeyError, json.JSONDecodeError) as exc:
            logger.warning("Skipping %s due to parse error: %s", job_dir.name, exc)

    if not rows:
        return _empty_dataframe()

    df = pd.DataFrame(rows, columns=_COLUMNS)
    return df


# Column order guaranteed by the output contract
_COLUMNS = [
    "job_name",
    "iptm",
    "ptm",
    "ranking_score",
    "fraction_disordered",
    "has_clash",
    "chain_iptm_mean",
    "cross_chain_pae_min",
    "nucleotide",
]


def _empty_dataframe() -> pd.DataFrame:
    return pd.DataFrame(columns=_COLUMNS)


if __name__ == "__main__":
    import sys

    logging.basicConfig(level=logging.INFO, stream=sys.stderr)
    df = parse_af3()
    if df.empty:
        print("No AF3 jobs found.")
    else:
        print(df.to_string())
        results_dir = Path(__file__).resolve().parent.parent.parent / "results"
        results_dir.mkdir(exist_ok=True)
        out_csv = results_dir / "af3_scores.csv"
        df.to_csv(out_csv, index=False)
        print(f"\nSaved to {out_csv}")
