"""
docking_plots.py — Matplotlib/Seaborn visualizations for virtual screening results.

Consumes DataFrames produced by :func:`~src.docking.triage.triage_hits` and
generates publication-quality PNGs to ``results/plots/docking/``.

All public functions return a ``matplotlib.figure.Figure`` for inline Jupyter
display and accept a *save* keyword to control whether the PNG is written.

Install dependencies:
    pip install matplotlib seaborn
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Optional matplotlib / seaborn imports
# ---------------------------------------------------------------------------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    _MPL_AVAILABLE = True
except ImportError:
    _MPL_AVAILABLE = False
    logger.warning(
        "matplotlib is not installed. Install it with:\n"
        "    pip install matplotlib\n"
        "All docking plot functions will return None until it is available."
    )

try:
    import seaborn as sns
    _SNS_AVAILABLE = True
except ImportError:
    _SNS_AVAILABLE = False
    logger.debug("seaborn not installed; basic matplotlib styles will be used.")

try:
    import numpy as np
    _NP_AVAILABLE = True
except ImportError:
    _NP_AVAILABLE = False
    logger.warning("numpy not installed; some curve-fitting features will be unavailable.")

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _docking_plot_dir() -> Path:
    """Return ``results/plots/docking/`` relative to the repo root, creating it."""
    repo_root = Path(__file__).resolve().parent.parent.parent
    d = repo_root / "results" / "plots" / "docking"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _save_fig(fig: "plt.Figure", filename: str) -> None:
    """Save *fig* as a PNG to the docking plot directory."""
    out = _docking_plot_dir() / filename
    fig.savefig(str(out), dpi=150, bbox_inches="tight")
    logger.info("Saved docking plot to %s", out)


def _apply_style() -> None:
    """Apply a clean seaborn style if available, else use matplotlib defaults."""
    if _SNS_AVAILABLE:
        sns.set_theme(style="whitegrid", font_scale=1.0)


def _no_hits_fig(message: str = "No hits found") -> Optional["plt.Figure"]:
    """Return a figure containing only a centred text message."""
    if not _MPL_AVAILABLE:
        return None
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.text(
        0.5, 0.5, message,
        ha="center", va="center",
        fontsize=14, transform=ax.transAxes,
    )
    ax.axis("off")
    fig.tight_layout()
    return fig


def _gnina_available(hits_df: pd.DataFrame) -> bool:
    """Return True if hits_df has a gnina_score column with at least one finite value."""
    if "gnina_score" not in hits_df.columns:
        return False
    return not hits_df["gnina_score"].isna().all()


# ---------------------------------------------------------------------------
# score_distribution
# ---------------------------------------------------------------------------

def plot_score_distribution(
    hits_df: pd.DataFrame,
    decoy_mean: float,
    decoy_std: float,
    save: bool = True,
) -> Optional["plt.Figure"]:
    """Histogram of Vina scores with a decoy normal-distribution overlay.

    Parameters
    ----------
    hits_df : pd.DataFrame
        DataFrame as returned by :func:`~src.docking.triage.triage_hits`.
        Must contain ``vina_score``.
    decoy_mean : float
        Mean of the decoy score distribution.
    decoy_std : float
        Standard deviation of the decoy score distribution.
    save : bool
        If True, write PNG to ``results/plots/docking/score_distribution.png``.

    Returns
    -------
    matplotlib.figure.Figure or None
    """
    if not _MPL_AVAILABLE:
        return None

    _apply_style()

    if hits_df.empty:
        fig = _no_hits_fig("No hits found")
        if save:
            _save_fig(fig, "score_distribution.png")
        return fig

    fig, ax = plt.subplots(figsize=(8, 5))

    scores = hits_df["vina_score"].dropna().values

    # Histogram of all vina scores (grey bars)
    n_bins = max(15, len(scores) // 5)
    ax.hist(scores, bins=n_bins, color="grey", alpha=0.7, label="Vina scores")

    # Overlay normal distribution curve from decoy_mean / decoy_std
    if _NP_AVAILABLE and decoy_std > 0:
        x_min = min(scores.min(), decoy_mean - 4 * decoy_std)
        x_max = max(scores.max(), decoy_mean + 4 * decoy_std)
        x = np.linspace(x_min, x_max, 300)
        pdf = (
            np.exp(-0.5 * ((x - decoy_mean) / decoy_std) ** 2)
            / (decoy_std * np.sqrt(2 * np.pi))
        )
        # Scale PDF to match histogram counts (density → counts)
        bin_width = (scores.max() - scores.min()) / n_bins
        ax.plot(
            x,
            pdf * len(scores) * bin_width,
            color="black",
            linestyle="--",
            linewidth=1.8,
            label=f"Decoy distribution\n(μ={decoy_mean:.2f}, σ={decoy_std:.2f})",
        )

    # Vertical dashed red line at Z = -2.5 cutoff
    cutoff_score = decoy_mean + (-2.5 * decoy_std)
    ax.axvline(
        x=cutoff_score,
        color="red",
        linestyle="--",
        linewidth=1.8,
        label=f"Z = -2.5 cutoff ({cutoff_score:.2f} kcal/mol)",
    )

    # Annotate number of hits below cutoff
    n_below = int((scores < cutoff_score).sum())
    ax.text(
        cutoff_score - 0.05,
        ax.get_ylim()[1] * 0.85 if ax.get_ylim()[1] > 0 else 1,
        f"n = {n_below} hits\nbelow cutoff",
        ha="right",
        va="top",
        fontsize=9,
        color="red",
    )

    ax.set_xlabel("Vina score (kcal/mol)", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title("Vina Score Distribution", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    fig.tight_layout()

    if save:
        _save_fig(fig, "score_distribution.png")
    return fig


# ---------------------------------------------------------------------------
# zscore_ranked
# ---------------------------------------------------------------------------

def plot_zscore_ranked(
    hits_df: pd.DataFrame,
    cutoff: float = -2.5,
    top_n: int = 50,
    save: bool = True,
) -> Optional["plt.Figure"]:
    """Bar chart of top-N compounds ranked by Vina Z-score.

    Parameters
    ----------
    hits_df : pd.DataFrame
        DataFrame from :func:`~src.docking.triage.triage_hits`.
        Must contain ``compound_id``, ``vina_zscore``, ``consensus_flag``.
    cutoff : float
        Z-score threshold shown as a horizontal dashed red line.
    top_n : int
        Number of top compounds (most negative Z-score) to display.
    save : bool
        If True, write PNG to ``results/plots/docking/zscore_ranked.png``.

    Returns
    -------
    matplotlib.figure.Figure or None
    """
    if not _MPL_AVAILABLE:
        return None

    _apply_style()

    if hits_df.empty:
        fig = _no_hits_fig("No hits found")
        if save:
            _save_fig(fig, "zscore_ranked.png")
        return fig

    # Sort by vina_zscore ascending (most negative first) and take top_n
    subset = (
        hits_df.dropna(subset=["vina_zscore"])
        .sort_values("vina_zscore", ascending=True)
        .head(top_n)
    )

    if subset.empty:
        fig = _no_hits_fig("No Z-score data available")
        if save:
            _save_fig(fig, "zscore_ranked.png")
        return fig

    bar_colors = [
        "green" if flag else "grey"
        for flag in subset["consensus_flag"].tolist()
    ]

    fig_width = max(10, len(subset) * 0.4)
    fig, ax = plt.subplots(figsize=(fig_width, 5))

    ax.bar(
        range(len(subset)),
        subset["vina_zscore"].values,
        color=bar_colors,
        edgecolor="white",
        linewidth=0.5,
    )

    # Horizontal dashed red line at cutoff
    ax.axhline(y=cutoff, color="red", linestyle="--", linewidth=1.6,
               label=f"Cutoff (Z = {cutoff})")

    # Legend patches
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="green", label="Consensus hit"),
        Patch(facecolor="grey", label="Non-consensus"),
        plt.Line2D([0], [0], color="red", linestyle="--", linewidth=1.6,
                   label=f"Cutoff (Z = {cutoff})"),
    ]
    ax.legend(handles=legend_elements, fontsize=9)

    ax.set_xticks(range(len(subset)))
    ax.set_xticklabels(subset["compound_id"].tolist(), rotation=90, fontsize=8)
    ax.set_xlabel("Compound ID", fontsize=11)
    ax.set_ylabel("Vina Z-score", fontsize=11)
    ax.set_title(
        f"Top {len(subset)} Compounds by Vina Z-score", fontsize=13, fontweight="bold"
    )
    fig.tight_layout()

    if save:
        _save_fig(fig, "zscore_ranked.png")
    return fig


# ---------------------------------------------------------------------------
# vina_gnina_scatter
# ---------------------------------------------------------------------------

def plot_vina_gnina_scatter(
    hits_df: pd.DataFrame,
    cutoff: float = -2.5,
    save: bool = True,
) -> Optional["plt.Figure"]:
    """Scatter plot of Vina Z-score vs GNINA Z-score.

    Parameters
    ----------
    hits_df : pd.DataFrame
        DataFrame from :func:`~src.docking.triage.triage_hits`.
        Must contain ``vina_zscore``, ``gnina_zscore``, ``consensus_flag``,
        ``compound_id``.
    cutoff : float
        Dashed quadrant lines drawn at this Z-score on both axes.
    save : bool
        If True, write PNG to ``results/plots/docking/vina_gnina_scatter.png``.

    Returns
    -------
    matplotlib.figure.Figure or None
    """
    if not _MPL_AVAILABLE:
        return None

    _apply_style()

    if hits_df.empty:
        fig = _no_hits_fig("No hits found")
        if save:
            _save_fig(fig, "vina_gnina_scatter.png")
        return fig

    # If GNINA scores are all None: show informational message plot
    if not _gnina_available(hits_df):
        fig = _no_hits_fig("GNINA scores not available")
        if save:
            _save_fig(fig, "vina_gnina_scatter.png")
        return fig

    subset = hits_df.dropna(subset=["vina_zscore", "gnina_zscore"]).copy()

    if subset.empty:
        fig = _no_hits_fig("No paired Vina/GNINA Z-score data")
        if save:
            _save_fig(fig, "vina_gnina_scatter.png")
        return fig

    fig, ax = plt.subplots(figsize=(8, 7))

    consensus = subset[subset["consensus_flag"] == True]
    non_consensus = subset[subset["consensus_flag"] != True]

    if not non_consensus.empty:
        ax.scatter(
            non_consensus["vina_zscore"],
            non_consensus["gnina_zscore"],
            color="grey",
            alpha=0.6,
            s=40,
            label="Non-consensus",
            zorder=2,
        )

    if not consensus.empty:
        ax.scatter(
            consensus["vina_zscore"],
            consensus["gnina_zscore"],
            color="green",
            alpha=0.8,
            s=60,
            label="Consensus hit",
            zorder=3,
        )

    # Quadrant lines
    ax.axvline(x=cutoff, color="red", linestyle="--", linewidth=1.4,
               label=f"Cutoff (Z = {cutoff})")
    ax.axhline(y=cutoff, color="red", linestyle="--", linewidth=1.4)

    # Annotate top 10 consensus hits with compound_id
    top_consensus = (
        consensus.sort_values("vina_zscore", ascending=True).head(10)
        if not consensus.empty else pd.DataFrame()
    )
    for _, row in top_consensus.iterrows():
        ax.annotate(
            str(row["compound_id"]),
            xy=(row["vina_zscore"], row["gnina_zscore"]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7,
            color="darkgreen",
        )

    ax.set_xlabel("Vina Z-score", fontsize=12)
    ax.set_ylabel("GNINA Z-score", fontsize=12)
    ax.set_title("Vina vs GNINA Z-score", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    fig.tight_layout()

    if save:
        _save_fig(fig, "vina_gnina_scatter.png")
    return fig


# ---------------------------------------------------------------------------
# top_hits
# ---------------------------------------------------------------------------

def plot_top_hits(
    hits_df: pd.DataFrame,
    top_n: int = 20,
    save: bool = True,
) -> Optional["plt.Figure"]:
    """Horizontal bar chart of top-N hits with Vina and GNINA Z-scores.

    Parameters
    ----------
    hits_df : pd.DataFrame
        DataFrame from :func:`~src.docking.triage.triage_hits`.
    top_n : int
        Number of top compounds to display.
    save : bool
        If True, write PNG to ``results/plots/docking/top_hits.png``.

    Returns
    -------
    matplotlib.figure.Figure or None
    """
    if not _MPL_AVAILABLE:
        return None

    _apply_style()

    if hits_df.empty:
        fig = _no_hits_fig("No hits found")
        if save:
            _save_fig(fig, "top_hits.png")
        return fig

    subset = (
        hits_df.dropna(subset=["vina_zscore"])
        .sort_values("vina_zscore", ascending=True)
        .head(top_n)
    )

    if subset.empty:
        fig = _no_hits_fig("No Z-score data available")
        if save:
            _save_fig(fig, "top_hits.png")
        return fig

    has_gnina = _gnina_available(subset)
    n = len(subset)

    # Build y-axis labels; prepend star for consensus hits
    labels = [
        f"★ {cid}" if flag else str(cid)
        for cid, flag in zip(subset["compound_id"], subset["consensus_flag"])
    ]

    fig_height = max(5, n * 0.45)
    fig, ax = plt.subplots(figsize=(10, fig_height))

    y_positions = list(range(n))

    if has_gnina:
        bar_height = 0.35
        vina_positions = [y + bar_height / 2 for y in y_positions]
        gnina_positions = [y - bar_height / 2 for y in y_positions]

        gnina_zscores = (
            subset["gnina_zscore"]
            .fillna(0)
            .values
        )

        bars_vina = ax.barh(
            vina_positions,
            subset["vina_zscore"].values,
            height=bar_height,
            color="#4472C4",
            label="Vina Z-score",
        )
        bars_gnina = ax.barh(
            gnina_positions,
            gnina_zscores,
            height=bar_height,
            color="#ED7D31",
            label="GNINA Z-score",
        )

        # Annotate bars with raw kcal/mol scores
        for bar, score in zip(bars_vina, subset["vina_score"].values):
            width = bar.get_width()
            ax.text(
                width - 0.05 if width < 0 else width + 0.05,
                bar.get_y() + bar.get_height() / 2,
                f"{score:.2f}",
                ha="right" if width < 0 else "left",
                va="center",
                fontsize=7,
                color="#4472C4",
            )

        gnina_raw = subset["gnina_score"].values
        for bar, score in zip(bars_gnina, gnina_raw):
            if score is None or (hasattr(score, "__class__") and
                                  score.__class__.__name__ == "float" and
                                  _NP_AVAILABLE and np.isnan(score)):
                continue
            width = bar.get_width()
            ax.text(
                width - 0.05 if width < 0 else width + 0.05,
                bar.get_y() + bar.get_height() / 2,
                f"{float(score):.2f}",
                ha="right" if width < 0 else "left",
                va="center",
                fontsize=7,
                color="#ED7D31",
            )

        ax.set_yticks(y_positions)
        ax.set_yticklabels(labels, fontsize=9)
        ax.legend(fontsize=9)

    else:
        # Single bar (Vina only)
        bars = ax.barh(
            y_positions,
            subset["vina_zscore"].values,
            color="#4472C4",
            label="Vina Z-score",
        )
        for bar, score in zip(bars, subset["vina_score"].values):
            width = bar.get_width()
            ax.text(
                width - 0.05 if width < 0 else width + 0.05,
                bar.get_y() + bar.get_height() / 2,
                f"{score:.2f}",
                ha="right" if width < 0 else "left",
                va="center",
                fontsize=8,
                color="#4472C4",
            )
        ax.set_yticks(y_positions)
        ax.set_yticklabels(labels, fontsize=9)
        ax.legend(fontsize=9)

    ax.set_xlabel("Z-score", fontsize=11)
    ax.set_title(f"Top {n} Docking Hits", fontsize=13, fontweight="bold")
    ax.invert_yaxis()  # best score at top
    fig.tight_layout()

    if save:
        _save_fig(fig, "top_hits.png")
    return fig


# ---------------------------------------------------------------------------
# screening_summary (2x2 composite)
# ---------------------------------------------------------------------------

def plot_screening_summary(
    hits_df: pd.DataFrame,
    decoy_mean: float,
    decoy_std: float,
    save: bool = True,
) -> Optional["plt.Figure"]:
    """2x2 composite overview of the full virtual screen.

    Subplots::

        [0,0] score_distribution
        [0,1] zscore_ranked (top 30)
        [1,0] vina_gnina_scatter
        [1,1] top_hits (top 10)

    Parameters
    ----------
    hits_df : pd.DataFrame
        DataFrame from :func:`~src.docking.triage.triage_hits`.
    decoy_mean : float
        Mean decoy Vina score.
    decoy_std : float
        Standard deviation of decoy Vina scores.
    save : bool
        If True, write PNG to ``results/plots/docking/screening_summary.png``.

    Returns
    -------
    matplotlib.figure.Figure or None
    """
    if not _MPL_AVAILABLE:
        return None

    _apply_style()

    if hits_df.empty:
        fig = _no_hits_fig("No hits found")
        if save:
            _save_fig(fig, "screening_summary.png")
        return fig

    # Generate the four component figures (without saving)
    fig_sd = plot_score_distribution(hits_df, decoy_mean, decoy_std, save=False)
    fig_zr = plot_zscore_ranked(hits_df, top_n=30, save=False)
    fig_sc = plot_vina_gnina_scatter(hits_df, save=False)
    fig_th = plot_top_hits(hits_df, top_n=10, save=False)

    component_figs = [fig_sd, fig_zr, fig_sc, fig_th]
    titles = [
        "Score Distribution",
        "Z-score Ranked (top 30)",
        "Vina vs GNINA Z-score",
        "Top 10 Hits",
    ]

    fig, axes = plt.subplots(2, 2, figsize=(18, 12))
    axes_flat = axes.flatten()

    for idx, (comp_fig, title) in enumerate(zip(component_figs, titles)):
        ax_target = axes_flat[idx]

        if comp_fig is None:
            ax_target.text(
                0.5, 0.5, f"{title}\n(unavailable)",
                ha="center", va="center", transform=ax_target.transAxes,
                fontsize=10,
            )
            ax_target.axis("off")
            continue

        # Transfer artists from the component figure's first axes into the
        # summary axes by re-rendering into a raster image and embedding it.
        try:
            import io
            buf = io.BytesIO()
            comp_fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
            buf.seek(0)
            img = plt.imread(buf)
            ax_target.imshow(img, aspect="auto")
            ax_target.axis("off")
            ax_target.set_title(title, fontsize=11, fontweight="bold", pad=6)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Could not embed subplot %s: %s", title, exc)
            ax_target.text(
                0.5, 0.5, f"{title}\n(render error)",
                ha="center", va="center", transform=ax_target.transAxes,
                fontsize=10,
            )
            ax_target.axis("off")
        finally:
            plt.close(comp_fig)

    fig.suptitle("Virtual Screening Summary", fontsize=15, fontweight="bold", y=1.01)
    fig.tight_layout()

    if save:
        _save_fig(fig, "screening_summary.png")
    return fig


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    print("=== docking_plots.py self-test ===")

    # Build a synthetic triage DataFrame
    _rng = None
    if _NP_AVAILABLE:
        _rng = np.random.default_rng(42)
        n = 40
        vina_scores = _rng.normal(-7.0, 1.5, n).tolist()
        gnina_scores = [v + _rng.normal(0, 0.5) for v in vina_scores]
        decoy_m, decoy_s = -6.0, 1.0

        from src.docking.triage import compute_zscores
        import pandas as _pd
        vz = compute_zscores(_pd.Series(vina_scores), decoy_m, decoy_s).tolist()
        gz = compute_zscores(_pd.Series(gnina_scores), decoy_m, decoy_s).tolist()

        df = _pd.DataFrame({
            "compound_id": [f"cpd_{i:03d}" for i in range(n)],
            "smiles": ["CCO"] * n,
            "vina_score": vina_scores,
            "gnina_score": gnina_scores,
            "vina_zscore": vz,
            "gnina_zscore": gz,
            "pose_file": [f"/tmp/cpd_{i:03d}_pose.pdbqt" for i in range(n)],
            "consensus_flag": [v <= -2.5 and g <= -2.5 for v, g in zip(vz, gz)],
        })

        print(f"Synthetic data: {len(df)} compounds, "
              f"{df['consensus_flag'].sum()} consensus hits")

        figs = {
            "score_distribution": plot_score_distribution(df, decoy_m, decoy_s),
            "zscore_ranked":      plot_zscore_ranked(df),
            "vina_gnina_scatter": plot_vina_gnina_scatter(df),
            "top_hits":           plot_top_hits(df),
            "screening_summary":  plot_screening_summary(df, decoy_m, decoy_s),
        }
        for name, fig in figs.items():
            status = "OK" if fig is not None else "FAILED"
            print(f"  {name}: {status}")
            if fig is not None:
                plt.close(fig)

        # Test GNINA-null path
        df_no_gnina = df.copy()
        df_no_gnina["gnina_score"] = None
        df_no_gnina["gnina_zscore"] = float("nan")
        fig_no_g = plot_vina_gnina_scatter(df_no_gnina)
        print(f"  vina_gnina_scatter (no GNINA): {'OK' if fig_no_g is not None else 'FAILED'}")
        if fig_no_g is not None:
            plt.close(fig_no_g)

        # Test empty DataFrame path
        empty_df = _pd.DataFrame(columns=df.columns)
        fig_empty = plot_score_distribution(empty_df, decoy_m, decoy_s)
        print(f"  score_distribution (empty df): {'OK' if fig_empty is not None else 'FAILED'}")
        if fig_empty is not None:
            plt.close(fig_empty)

    else:
        print("numpy not available — skipping data generation tests")

    print("=== self-test complete ===")
