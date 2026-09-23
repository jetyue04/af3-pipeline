"""AlphaFold 3 confidence score visualizations using matplotlib/seaborn.

All functions save PNGs to results/plots/af3/ and return fig objects so
they also display inline in Jupyter notebooks.

Install dependencies:
    pip install matplotlib seaborn pandas numpy
"""

import json
import logging
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Optional imports
# ---------------------------------------------------------------------------
try:
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib import cm
    from matplotlib.colors import LinearSegmentedColormap
    _MPL_AVAILABLE = True
except ImportError:  # pragma: no cover
    _MPL_AVAILABLE = False
    logger.warning(
        "matplotlib is not installed. All plot functions will return None. "
        "Install with: pip install matplotlib"
    )

try:
    import seaborn as sns
    _SNS_AVAILABLE = True
except ImportError:  # pragma: no cover
    _SNS_AVAILABLE = False
    logger.warning(
        "seaborn is not installed. Heatmaps will fall back to matplotlib imshow. "
        "Install with: pip install seaborn"
    )

# AF3 pLDDT colour scheme (matches AlphaFold colouring)
PLDDT_COLORS = {
    "very_high": "#0053D6",  # >= 90  dark blue
    "confident": "#65CBF3",  # 70-90  light blue
    "low":       "#FFDB13",  # 50-70  yellow
    "very_low":  "#FF7D45",  # < 50   orange
}

# Palette for GTP / GDP comparisons
COLOR_GTP = "#4C72B0"  # blue
COLOR_GDP = "#DD8452"  # orange


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _repo_root() -> Path:
    """Return the repository root (three levels above this file)."""
    return Path(__file__).resolve().parent.parent.parent


def _af3_plot_dir() -> Path:
    """Return results/plots/af3/, creating it if needed."""
    d = _repo_root() / "results" / "plots" / "af3"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _check_mpl(func_name: str) -> bool:
    """Return True if matplotlib is available, else log an error."""
    if not _MPL_AVAILABLE:
        logger.error(
            "%s requires matplotlib. Install with: pip install matplotlib",
            func_name,
        )
        return False
    return True


def _load_confidences(confidences_json_path: str) -> dict:
    """Load and return the full _confidences.json as a dict."""
    with open(confidences_json_path, "r") as fh:
        return json.load(fh)


def _get_chain_lengths(data: dict) -> list[int]:
    """Extract per-chain residue counts from confidences data.

    Tries ``token_chain_ids`` (list of chain letters, one per residue),
    then ``chain_lengths``, then falls back to the whole PAE matrix size.
    """
    # Preferred: token_chain_ids lists one entry per residue token
    if "token_chain_ids" in data:
        ids = data["token_chain_ids"]
        lengths: dict[str, int] = {}
        order: list[str] = []
        for c in ids:
            if c not in lengths:
                lengths[c] = 0
                order.append(c)
            lengths[c] += 1
        return [lengths[c] for c in order]

    # Second choice: explicit chain_lengths field
    if "chain_lengths" in data:
        return list(data["chain_lengths"])

    # Last resort: derive from PAE matrix dimension
    pae = data.get("pae") or data.get("predicted_aligned_error")
    if pae is not None:
        n = len(pae)
        return [n]  # treat as single chain of length n

    return []


def _chain_boundaries(chain_lengths: list[int]) -> list[int]:
    """Return cumulative residue indices marking chain end boundaries."""
    boundaries: list[int] = []
    cumsum = 0
    for length in chain_lengths[:-1]:  # no line after last chain
        cumsum += length
        boundaries.append(cumsum)
    return boundaries


def _get_plddt(data: dict) -> list[float]:
    """Extract per-residue pLDDT from confidences data.

    Tries ``atom_plddts`` first (AF3 ≥ v2), then ``plddt`` (AF3 v1 / summary).
    Averages atom-level values to residue-level using ``atom_chain_ids`` /
    ``token_chain_ids`` when available; otherwise returns raw values.
    """
    # atom_plddts: one value per heavy atom; we need per-residue averages
    if "atom_plddts" in data:
        atom_plddts = np.array(data["atom_plddts"], dtype=float)

        # Try to group atoms to residues via atom_res_id or token_res_ids
        if "atom_res_id" in data:
            res_ids = np.array(data["atom_res_id"])
            unique_ids = []
            seen = set()
            for r in res_ids:
                if r not in seen:
                    unique_ids.append(r)
                    seen.add(r)
            residue_plddt = [
                float(atom_plddts[res_ids == r].mean()) for r in unique_ids
            ]
            return residue_plddt

        # Fallback: return atom-level values as-is (still valid for plotting)
        return list(atom_plddts)

    if "plddt" in data:
        return list(data["plddt"])

    return []


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def plot_pae_heatmap(
    confidences_json_path: str,
    job_name: str,
    save: bool = True,
) -> "plt.Figure | None":
    """Plot the PAE matrix from a _confidences.json file as a seaborn heatmap.

    Parameters
    ----------
    confidences_json_path:
        Path to the AF3 full ``*_confidences.json`` file (not the summary).
    job_name:
        Job name used in the plot title and output filename.
    save:
        If True, save the figure to results/plots/af3/<job_name>_pae_heatmap.png.

    Returns
    -------
    matplotlib.figure.Figure or None (if matplotlib is unavailable).
    """
    if not _check_mpl("plot_pae_heatmap"):
        return None

    data = _load_confidences(confidences_json_path)

    # Accept both 'pae' and 'predicted_aligned_error' field names
    pae_raw = data.get("pae") or data.get("predicted_aligned_error")
    if pae_raw is None:
        logger.error(
            "plot_pae_heatmap: 'pae' / 'predicted_aligned_error' key not found in %s",
            confidences_json_path,
        )
        return None

    pae_matrix = np.array(pae_raw, dtype=float)
    n = pae_matrix.shape[0]

    chain_lengths = _get_chain_lengths(data)
    boundaries = _chain_boundaries(chain_lengths)

    # Build colour map: dark blue (0) -> white (30)
    pae_cmap = LinearSegmentedColormap.from_list(
        "pae_cmap",
        [(0.0, "#0A0A6B"), (0.5, "#5BA3D9"), (1.0, "#FFFFFF")],
    )

    fig, ax = plt.subplots(figsize=(8, 7))

    if _SNS_AVAILABLE:
        sns.heatmap(
            pae_matrix,
            ax=ax,
            cmap=pae_cmap,
            vmin=0,
            vmax=30,
            xticklabels=False,
            yticklabels=False,
            cbar_kws={"label": "PAE (Å)"},
            square=True,
        )
    else:
        im = ax.imshow(pae_matrix, cmap=pae_cmap, vmin=0, vmax=30, aspect="auto")
        fig.colorbar(im, ax=ax, label="PAE (Å)")
        ax.set_xticks([])
        ax.set_yticks([])

    # Chain boundary lines
    for b in boundaries:
        ax.axvline(x=b, color="black", linestyle="--", linewidth=1.0, alpha=0.7)
        ax.axhline(y=b, color="black", linestyle="--", linewidth=1.0, alpha=0.7)

    ax.set_title(f"PAE — {job_name}", fontsize=13, fontweight="bold")
    ax.set_xlabel("Aligned residue")
    ax.set_ylabel("Scored residue")

    fig.tight_layout()

    if save:
        out_path = _af3_plot_dir() / f"{job_name}_pae_heatmap.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        logger.info("Saved PAE heatmap: %s", out_path)

    return fig


def plot_plddt_per_residue(
    confidences_json_path: str,
    job_name: str,
    save: bool = True,
) -> "plt.Figure | None":
    """Plot per-residue pLDDT as a colour-coded line chart.

    Parameters
    ----------
    confidences_json_path:
        Path to the AF3 full ``*_confidences.json`` file.
    job_name:
        Job name used in the plot title and output filename.
    save:
        If True, save to results/plots/af3/<job_name>_plddt.png.

    Returns
    -------
    matplotlib.figure.Figure or None.
    """
    if not _check_mpl("plot_plddt_per_residue"):
        return None

    data = _load_confidences(confidences_json_path)
    plddt = _get_plddt(data)

    if not plddt:
        logger.error(
            "plot_plddt_per_residue: no pLDDT data found in %s",
            confidences_json_path,
        )
        return None

    plddt = np.array(plddt, dtype=float)
    n = len(plddt)
    residues = np.arange(n)

    chain_lengths = _get_chain_lengths(data)
    boundaries = _chain_boundaries(chain_lengths)

    fig, ax = plt.subplots(figsize=(max(8, n // 40), 4))

    # Shaded confidence bands (drawn first, behind the line)
    ax.axhspan(90, 100, color=PLDDT_COLORS["very_high"], alpha=0.15, zorder=0)
    ax.axhspan(70,  90, color=PLDDT_COLORS["confident"],  alpha=0.15, zorder=0)
    ax.axhspan(50,  70, color=PLDDT_COLORS["low"],         alpha=0.15, zorder=0)
    ax.axhspan(0,   50, color=PLDDT_COLORS["very_low"],    alpha=0.15, zorder=0)

    # Colour-segment the line by pLDDT band
    def _band_color(v: float) -> str:
        if v >= 90:
            return PLDDT_COLORS["very_high"]
        if v >= 70:
            return PLDDT_COLORS["confident"]
        if v >= 50:
            return PLDDT_COLORS["low"]
        return PLDDT_COLORS["very_low"]

    # Plot one segment per contiguous same-band run
    start = 0
    current_color = _band_color(plddt[0])
    for i in range(1, n):
        c = _band_color(plddt[i])
        if c != current_color or i == n - 1:
            end = i + 1 if i == n - 1 else i
            ax.plot(
                residues[start:end],
                plddt[start:end],
                color=current_color,
                linewidth=1.5,
                zorder=2,
            )
            start = i
            current_color = c

    # Horizontal reference lines
    for threshold, label in [(90, "Very high (≥90)"), (70, "Confident (≥70)"), (50, "Low (≥50)")]:
        ax.axhline(
            y=threshold,
            color="grey",
            linestyle="--",
            linewidth=0.8,
            alpha=0.6,
            zorder=1,
        )
        ax.text(
            n * 0.01, threshold + 1, label,
            fontsize=7, color="grey", va="bottom",
        )

    # Chain boundary lines (multi-chain)
    for b in boundaries:
        ax.axvline(x=b, color="black", linestyle="--", linewidth=1.0, alpha=0.6)

    ax.set_xlim(0, n - 1)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Residue index")
    ax.set_ylabel("pLDDT")
    ax.set_title(f"pLDDT per residue — {job_name}", fontsize=13, fontweight="bold")

    # Legend patches
    legend_patches = [
        mpatches.Patch(color=PLDDT_COLORS["very_high"], label="Very high ≥90"),
        mpatches.Patch(color=PLDDT_COLORS["confident"],  label="Confident 70-90"),
        mpatches.Patch(color=PLDDT_COLORS["low"],         label="Low 50-70"),
        mpatches.Patch(color=PLDDT_COLORS["very_low"],    label="Very low <50"),
    ]
    ax.legend(handles=legend_patches, loc="lower right", fontsize=8, framealpha=0.8)

    fig.tight_layout()

    if save:
        out_path = _af3_plot_dir() / f"{job_name}_plddt.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        logger.info("Saved pLDDT plot: %s", out_path)

    return fig


def plot_dg_comparison(
    scores_df,
    pairs: list[tuple[str, str]],
    save: bool = True,
) -> "plt.Figure | None":
    """Side-by-side iptm bar chart comparing GTP vs GDP states for given pairs.

    Parameters
    ----------
    scores_df:
        DataFrame returned by ``parse_af3()``.  Required columns:
        ``job_name``, ``iptm``.
    pairs:
        List of ``(gtp_job_name, gdp_job_name)`` tuples.
    save:
        If True, save to results/plots/af3/dg_comparison.png.

    Returns
    -------
    matplotlib.figure.Figure or None.
    """
    if not _check_mpl("plot_dg_comparison"):
        return None

    if not pairs:
        logger.warning("plot_dg_comparison: pairs list is empty.")
        return None

    n_pairs = len(pairs)
    fig, ax = plt.subplots(figsize=(max(6, n_pairs * 2.5), 5))

    bar_width = 0.35
    x = np.arange(n_pairs)

    # Index scores_df by job_name for fast lookup
    score_lookup = {}
    if scores_df is not None and not scores_df.empty:
        score_lookup = scores_df.set_index("job_name")["iptm"].to_dict()

    gtp_iptm = [score_lookup.get(g, float("nan")) for g, _ in pairs]
    gdp_iptm = [score_lookup.get(d, float("nan")) for _, d in pairs]

    bars_gtp = ax.bar(
        x - bar_width / 2, gtp_iptm,
        width=bar_width, color=COLOR_GTP, label="GTP", zorder=3,
    )
    bars_gdp = ax.bar(
        x + bar_width / 2, gdp_iptm,
        width=bar_width, color=COLOR_GDP, label="GDP", zorder=3,
    )

    # Annotate bars
    for bar in list(bars_gtp) + list(bars_gdp):
        height = bar.get_height()
        if not np.isnan(height):
            ax.annotate(
                f"{height:.3f}",
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center", va="bottom",
                fontsize=8,
            )

    # Minimum credible interface threshold line
    ax.axhline(
        y=0.5,
        color="red",
        linestyle="--",
        linewidth=1.2,
        alpha=0.8,
        label="iptm = 0.5 (min credible)",
        zorder=2,
    )

    # x-tick labels: use the common prefix of the pair (strip nucleotide suffix)
    pair_labels = []
    for gtp_name, gdp_name in pairs:
        # Use GTP name without trailing nucleotide info as label
        label = gtp_name.replace("_gtp", "").replace("_GTP", "")
        pair_labels.append(label)

    ax.set_xticks(x)
    ax.set_xticklabels(pair_labels, rotation=20, ha="right", fontsize=9)
    ax.set_ylabel("ipTM")
    ax.set_ylim(0, 1.05)
    ax.set_title("GTP vs GDP ipTM comparison", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(axis="y", alpha=0.4, zorder=0)

    fig.tight_layout()

    if save:
        out_path = _af3_plot_dir() / "dg_comparison.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        logger.info("Saved DG comparison plot: %s", out_path)

    return fig


def plot_bsa_comparison(
    interface_results: dict,
    save: bool = True,
) -> "plt.Figure | None":
    """Grouped bar chart comparing BSA and contact count for GTP vs GDP pairs.

    Parameters
    ----------
    interface_results:
        ``{job_name: {"bsa": float, "contact_count": int, ...}}``
        as returned by ``compute_interface()``.
    save:
        If True, save to results/plots/af3/bsa_comparison.png.

    Returns
    -------
    matplotlib.figure.Figure or None.
    """
    if not _check_mpl("plot_bsa_comparison"):
        return None

    if not interface_results:
        logger.warning("plot_bsa_comparison: interface_results is empty.")
        return None

    job_names = list(interface_results.keys())

    # Separate GTP and GDP jobs, build groups by shared prefix
    gtp_jobs = {k: v for k, v in interface_results.items()
                if "gtp" in k.lower()}
    gdp_jobs = {k: v for k, v in interface_results.items()
                if "gdp" in k.lower()}
    other_jobs = {k: v for k, v in interface_results.items()
                  if "gtp" not in k.lower() and "gdp" not in k.lower()}

    # Build groups: each GTP job paired with a GDP job sharing the same base
    groups: list[tuple[str, str | None, str | None]] = []  # (label, gtp_key, gdp_key)
    used_gtp: set[str] = set()
    used_gdp: set[str] = set()

    for gk in gtp_jobs:
        base = gk.lower().replace("_gtp", "").replace("gtp_", "")
        matched_dk = None
        for dk in gdp_jobs:
            if dk not in used_gdp:
                dbase = dk.lower().replace("_gdp", "").replace("gdp_", "")
                if base == dbase:
                    matched_dk = dk
                    break
        groups.append((gk.replace("_gtp", "").replace("_GTP", ""), gk, matched_dk))
        used_gtp.add(gk)
        if matched_dk:
            used_gdp.add(matched_dk)

    # Add unmatched GDP and other jobs
    for dk in gdp_jobs:
        if dk not in used_gdp:
            groups.append((dk.replace("_gdp", "").replace("_GDP", ""), None, dk))
    for ok in other_jobs:
        groups.append((ok, ok, None))

    n_groups = len(groups)
    if n_groups == 0:
        logger.warning("plot_bsa_comparison: no groups could be formed.")
        return None

    bar_width = 0.35
    x = np.arange(n_groups)

    fig, ax1 = plt.subplots(figsize=(max(6, n_groups * 2.5), 5))
    ax2 = ax1.twinx()

    gtp_bsa, gdp_bsa = [], []
    gtp_cc, gdp_cc = [], []
    labels = []

    for label, gk, dk in groups:
        labels.append(label)
        gtp_bsa.append(interface_results[gk]["bsa"] if gk and gk in interface_results else float("nan"))
        gdp_bsa.append(interface_results[dk]["bsa"] if dk and dk in interface_results else float("nan"))
        gtp_cc.append(interface_results[gk]["contact_count"] if gk and gk in interface_results else float("nan"))
        gdp_cc.append(interface_results[dk]["contact_count"] if dk and dk in interface_results else float("nan"))

    bars_gtp = ax1.bar(
        x - bar_width / 2, gtp_bsa,
        width=bar_width, color=COLOR_GTP, alpha=0.85, label="GTP BSA", zorder=3,
    )
    bars_gdp = ax1.bar(
        x + bar_width / 2, gdp_bsa,
        width=bar_width, color=COLOR_GDP, alpha=0.85, label="GDP BSA", zorder=3,
    )

    # Contact count on secondary y-axis (line plot)
    ax2.plot(
        x - bar_width / 2, gtp_cc,
        marker="o", color=COLOR_GTP, linestyle="--",
        linewidth=1.5, markersize=6, label="GTP contacts", zorder=4,
    )
    ax2.plot(
        x + bar_width / 2, gdp_cc,
        marker="s", color=COLOR_GDP, linestyle="--",
        linewidth=1.5, markersize=6, label="GDP contacts", zorder=4,
    )

    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=20, ha="right", fontsize=9)
    ax1.set_ylabel("BSA (Å²)", color="black")
    ax2.set_ylabel("Contact count", color="dimgrey")
    ax1.set_title("BSA and contact count: GTP vs GDP", fontsize=13, fontweight="bold")
    ax1.grid(axis="y", alpha=0.3, zorder=0)

    # Combined legend
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles1 + handles2, labels1 + labels2, fontsize=8, loc="upper right")

    fig.tight_layout()

    if save:
        out_path = _af3_plot_dir() / "bsa_comparison.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        logger.info("Saved BSA comparison plot: %s", out_path)

    return fig


def plot_score_summary_table(
    scores_df,
    save: bool = True,
) -> "plt.Figure | None":
    """Render scores_df as a colour-coded matplotlib table.

    Colour coding:
    - iptm >= 0.75  → green
    - iptm 0.5-0.75 → yellow
    - iptm < 0.5    → red
    - fraction_disordered > 0.5 → orange

    Parameters
    ----------
    scores_df:
        DataFrame from ``parse_af3()``.
    save:
        If True, save to results/plots/af3/score_summary.png.

    Returns
    -------
    matplotlib.figure.Figure or None.
    """
    if not _check_mpl("plot_score_summary_table"):
        return None

    if scores_df is None or scores_df.empty:
        logger.warning("plot_score_summary_table: scores_df is empty.")
        return None

    import pandas as pd  # local import – pandas is an existing dependency

    df = scores_df.copy()

    # Select and order display columns (show only the most informative subset)
    display_cols = [
        c for c in [
            "job_name", "iptm", "ptm", "ranking_score",
            "fraction_disordered", "has_clash",
            "chain_iptm_mean", "cross_chain_pae_min", "nucleotide",
        ]
        if c in df.columns
    ]
    df = df[display_cols].reset_index(drop=True)

    n_rows, n_cols = df.shape

    # Build cell colours: default white
    cell_colors = [["white"] * n_cols for _ in range(n_rows)]

    iptm_col = display_cols.index("iptm") if "iptm" in display_cols else None
    fdis_col = (
        display_cols.index("fraction_disordered")
        if "fraction_disordered" in display_cols
        else None
    )

    for row_idx in range(n_rows):
        if iptm_col is not None:
            iptm_val = df.iloc[row_idx, iptm_col]
            try:
                iptm_val = float(iptm_val)
                if iptm_val >= 0.75:
                    cell_colors[row_idx][iptm_col] = "#A8D5A2"   # green
                elif iptm_val >= 0.5:
                    cell_colors[row_idx][iptm_col] = "#FFF3A3"   # yellow
                else:
                    cell_colors[row_idx][iptm_col] = "#F4A4A4"   # red
            except (TypeError, ValueError):
                pass

        if fdis_col is not None:
            fdis_val = df.iloc[row_idx, fdis_col]
            try:
                if float(fdis_val) > 0.5:
                    cell_colors[row_idx][fdis_col] = "#FFCC88"   # orange
            except (TypeError, ValueError):
                pass

    # Format float columns to 3 decimal places
    def _fmt(v):
        try:
            return f"{float(v):.3f}"
        except (TypeError, ValueError):
            return str(v) if v is not None else ""

    table_data = []
    for row_idx in range(n_rows):
        row = []
        for col_idx, col in enumerate(display_cols):
            val = df.iloc[row_idx, col_idx]
            if col in ("iptm", "ptm", "ranking_score", "fraction_disordered",
                       "chain_iptm_mean", "cross_chain_pae_min"):
                row.append(_fmt(val))
            else:
                row.append(str(val) if val is not None else "")
        table_data.append(row)

    # Adjust figure height to row count
    fig_height = max(2, 0.35 * (n_rows + 1) + 1.0)
    fig, ax = plt.subplots(figsize=(max(14, n_cols * 1.8), fig_height))
    ax.axis("off")

    col_widths = [max(0.08, min(0.25, len(c) * 0.012)) for c in display_cols]
    col_widths_norm = [w / sum(col_widths) for w in col_widths]

    tbl = ax.table(
        cellText=table_data,
        colLabels=display_cols,
        cellColours=cell_colors,
        colWidths=col_widths_norm,
        loc="center",
        cellLoc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1.0, 1.4)

    # Style header row
    for col_idx, _ in enumerate(display_cols):
        tbl[(0, col_idx)].set_facecolor("#2C3E50")
        tbl[(0, col_idx)].set_text_props(color="white", fontweight="bold")

    ax.set_title(
        "AF3 Score Summary",
        fontsize=13, fontweight="bold", pad=10,
    )

    # Legend
    legend_patches = [
        mpatches.Patch(color="#A8D5A2", label="iptm ≥ 0.75"),
        mpatches.Patch(color="#FFF3A3", label="iptm 0.5–0.75"),
        mpatches.Patch(color="#F4A4A4", label="iptm < 0.5"),
        mpatches.Patch(color="#FFCC88", label="fraction_disordered > 0.5"),
    ]
    fig.legend(
        handles=legend_patches,
        loc="lower center",
        ncol=4,
        fontsize=8,
        framealpha=0.8,
        bbox_to_anchor=(0.5, 0.01),
    )

    fig.tight_layout(rect=[0, 0.06, 1, 1])

    if save:
        out_path = _af3_plot_dir() / "score_summary.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        logger.info("Saved score summary table: %s", out_path)

    return fig
