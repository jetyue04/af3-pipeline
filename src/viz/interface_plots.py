"""Interface analysis visualizations using matplotlib/seaborn.

All functions save PNGs to results/plots/interface/ and return fig objects
so they also display inline in Jupyter notebooks.

Install dependencies:
    pip install matplotlib seaborn biopython numpy
"""

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
    from matplotlib.patches import Rectangle
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

# Colours
COLOR_GTP = "#4C72B0"
COLOR_GDP = "#DD8452"
COLOR_CONTACT = "#C44E52"   # red for contact residue highlights
COLOR_CA = "#AAAAAA"         # grey for all-CA scatter


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _repo_root() -> Path:
    """Return the repository root (three levels above this file)."""
    return Path(__file__).resolve().parent.parent.parent


def _interface_plot_dir() -> Path:
    """Return results/plots/interface/, creating it if needed."""
    d = _repo_root() / "results" / "plots" / "interface"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _check_mpl(func_name: str) -> bool:
    if not _MPL_AVAILABLE:
        logger.error(
            "%s requires matplotlib. Install with: pip install matplotlib",
            func_name,
        )
        return False
    return True


def _load_structure(cif_path: str):
    """Parse a CIF file and return the first Biopython model."""
    from Bio.PDB import MMCIFParser  # type: ignore
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("struct", cif_path)
    return next(iter(structure))


def _get_chain(model, chain_id: str):
    """Return a Biopython chain object by ID."""
    try:
        return model[chain_id]
    except KeyError:
        logger.warning("Chain '%s' not found in structure.", chain_id)
        return None


def _ca_coords_for_chain(chain) -> tuple[np.ndarray, list[int]]:
    """Return (Cα coordinate array, residue_seq_ids) for a Biopython chain.

    Returns a float array of shape (N, 3) and a list of 1-indexed sequence IDs.
    """
    coords: list[list[float]] = []
    seq_ids: list[int] = []
    for residue in chain.get_residues():
        if "CA" in residue:
            ca = residue["CA"]
            coords.append(ca.get_vector().get_array().tolist())
            seq_ids.append(residue.get_id()[1])
    return np.array(coords, dtype=float), seq_ids


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def plot_contact_heatmap(
    cif_path: str,
    chain_a_id: str,
    chain_b_id: str,
    job_name: str,
    cutoff: float = 5.0,
    save: bool = True,
) -> "plt.Figure | None":
    """Residue-residue Cα distance heatmap for the interface region.

    Calls ``interface.get_contact_residues()`` to find contact residues, then
    computes a full Cα-Cα distance matrix for the interface patch (±10
    residues around the first and last contact on each chain).  Contact cells
    (distance ≤ cutoff) are overlaid with black squares.

    Parameters
    ----------
    cif_path:
        Path to the AF3 CIF model file.
    chain_a_id:
        Chain identifier for chain A.
    chain_b_id:
        Chain identifier for chain B.
    job_name:
        Used in the title and output filename.
    cutoff:
        Distance cutoff (Å) for contact detection (default 5.0).
    save:
        If True, save to results/plots/interface/<job_name>_contact_heatmap.png.

    Returns
    -------
    matplotlib.figure.Figure or None.
    """
    if not _check_mpl("plot_contact_heatmap"):
        return None

    # Get contact residues (0-indexed from interface.get_contact_residues)
    from src.scoring.interface import get_contact_residues  # type: ignore

    contact_a_0idx, contact_b_0idx = get_contact_residues(
        cif_path, chain_a_id, chain_b_id, cutoff=cutoff
    )

    if not contact_a_0idx and not contact_b_0idx:
        logger.warning(
            "plot_contact_heatmap: no contacts found for %s at cutoff %.1f Å",
            job_name, cutoff,
        )

    # Load structure and extract Cα coordinates
    model = _load_structure(cif_path)
    chain_a = _get_chain(model, chain_a_id)
    chain_b = _get_chain(model, chain_b_id)

    if chain_a is None or chain_b is None:
        logger.error("plot_contact_heatmap: missing chain(s) in %s", cif_path)
        return None

    coords_a, seq_ids_a = _ca_coords_for_chain(chain_a)
    coords_b, seq_ids_b = _ca_coords_for_chain(chain_b)

    n_res_a = len(seq_ids_a)
    n_res_b = len(seq_ids_b)

    if n_res_a == 0 or n_res_b == 0:
        logger.error(
            "plot_contact_heatmap: no Cα atoms found in chain(s) %s/%s",
            chain_a_id, chain_b_id,
        )
        return None

    # Build ±10-residue window around contacts on each chain (0-indexed)
    WINDOW = 10

    def _window_slice(contact_indices: list[int], n_res: int) -> slice:
        if not contact_indices:
            return slice(0, n_res)
        lo = max(0, min(contact_indices) - WINDOW)
        hi = min(n_res, max(contact_indices) + WINDOW + 1)
        return slice(lo, hi)

    slice_a = _window_slice(contact_a_0idx, n_res_a)
    slice_b = _window_slice(contact_b_0idx, n_res_b)

    patch_coords_a = coords_a[slice_a]
    patch_coords_b = coords_b[slice_b]
    patch_seq_a = seq_ids_a[slice_a.start:slice_a.stop]
    patch_seq_b = seq_ids_b[slice_b.start:slice_b.stop]

    n_a = len(patch_coords_a)
    n_b = len(patch_coords_b)

    if n_a == 0 or n_b == 0:
        logger.error("plot_contact_heatmap: interface patch is empty.")
        return None

    # Compute Cα-Cα distance matrix
    dist_matrix = np.sqrt(
        ((patch_coords_a[:, None, :] - patch_coords_b[None, :, :]) ** 2).sum(axis=2)
    )  # shape (n_a, n_b)

    # Build set of contact cells in patch-local coordinates
    contact_set: set[tuple[int, int]] = set()
    offset_a = slice_a.start
    offset_b = slice_b.start
    for ia in contact_a_0idx:
        for ib in contact_b_0idx:
            la = ia - offset_a
            lb = ib - offset_b
            if 0 <= la < n_a and 0 <= lb < n_b:
                if dist_matrix[la, lb] <= cutoff:
                    contact_set.add((la, lb))

    # ------ Plot ------
    fig, ax = plt.subplots(figsize=(max(6, n_b // 3 + 2), max(5, n_a // 3 + 2)))

    cmap_dist = "coolwarm_r"
    vmin = 0.0
    vmax = min(30.0, dist_matrix.max() * 1.1)

    if _SNS_AVAILABLE:
        sns.heatmap(
            dist_matrix,
            ax=ax,
            cmap=cmap_dist,
            vmin=vmin,
            vmax=vmax,
            xticklabels=False,
            yticklabels=False,
            cbar_kws={"label": "Cα–Cα distance (Å)"},
        )
    else:
        im = ax.imshow(
            dist_matrix, cmap=cmap_dist, vmin=vmin, vmax=vmax, aspect="auto"
        )
        fig.colorbar(im, ax=ax, label="Cα–Cα distance (Å)")
        ax.set_xticks([])
        ax.set_yticks([])

    # Overlay contact squares
    for (row_i, col_j) in contact_set:
        ax.add_patch(
            Rectangle(
                (col_j, row_i), 1, 1,
                linewidth=1.2,
                edgecolor="black",
                facecolor="none",
                zorder=3,
            )
        )

    # Tick labels every 5 residues
    step_a = max(1, n_a // 10)
    step_b = max(1, n_b // 10)
    ax.set_yticks(np.arange(0, n_a, step_a) + 0.5)
    ax.set_yticklabels(
        [str(patch_seq_a[i]) for i in range(0, n_a, step_a)],
        fontsize=7,
    )
    ax.set_xticks(np.arange(0, n_b, step_b) + 0.5)
    ax.set_xticklabels(
        [str(patch_seq_b[i]) for i in range(0, n_b, step_b)],
        rotation=45, ha="right", fontsize=7,
    )

    ax.set_xlabel(f"Chain {chain_b_id} residue")
    ax.set_ylabel(f"Chain {chain_a_id} residue")
    ax.set_title(
        f"Interface contact heatmap — {job_name}\n"
        f"(black squares: Cα ≤ {cutoff:.1f} Å)",
        fontsize=11, fontweight="bold",
    )

    fig.tight_layout()

    if save:
        out_path = _interface_plot_dir() / f"{job_name}_contact_heatmap.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        logger.info("Saved contact heatmap: %s", out_path)

    return fig


def plot_grid_box_2d(
    cif_path: str,
    chain_id: str,
    grid: dict,
    contact_residues: list[int],
    job_name: str,
    save: bool = True,
) -> "plt.Figure | None":
    """2D x-y projection of receptor with docking grid box overlay.

    Plots all Cα atoms as small grey dots, highlights contact residues as
    larger coloured dots, and draws the grid box as a rectangle in the x-y
    plane.

    Parameters
    ----------
    cif_path:
        Path to the AF3 CIF model file.
    chain_id:
        Chain ID of the receptor chain to project.
    grid:
        Dictionary with keys ``center_x``, ``center_y``, ``center_z``,
        ``size_x``, ``size_y``, ``size_z`` (all in Å).
    contact_residues:
        List of 0-indexed residue positions to highlight.
    job_name:
        Used in the title and output filename.
    save:
        If True, save to results/plots/interface/<job_name>_grid_box.png.

    Returns
    -------
    matplotlib.figure.Figure or None.
    """
    if not _check_mpl("plot_grid_box_2d"):
        return None

    model = _load_structure(cif_path)
    chain = _get_chain(model, chain_id)

    if chain is None:
        logger.error("plot_grid_box_2d: chain %s not found in %s", chain_id, cif_path)
        return None

    coords, seq_ids = _ca_coords_for_chain(chain)

    if len(coords) == 0:
        logger.error("plot_grid_box_2d: no Cα atoms found for chain %s", chain_id)
        return None

    all_x = coords[:, 0]
    all_y = coords[:, 1]

    # Convert 0-indexed contact residues to sequence IDs for lookup
    contact_seq_ids = set()
    for idx in contact_residues:
        if 0 <= idx < len(seq_ids):
            contact_seq_ids.add(seq_ids[idx])

    fig, ax = plt.subplots(figsize=(7, 6))

    # All Cα scatter
    ax.scatter(
        all_x, all_y,
        s=10, color=COLOR_CA, alpha=0.6, label="All Cα", zorder=2,
    )

    # Contact residues highlighted
    c_mask = np.array([sid in contact_seq_ids for sid in seq_ids])
    if c_mask.any():
        ax.scatter(
            all_x[c_mask], all_y[c_mask],
            s=60, color=COLOR_CONTACT, alpha=0.9,
            label=f"Contact residues (n={c_mask.sum()})",
            zorder=3,
        )

    # Grid box rectangle (x-y projection)
    cx = grid["center_x"]
    cy = grid["center_y"]
    hx = grid["size_x"] / 2.0
    hy = grid["size_y"] / 2.0

    box_rect = Rectangle(
        (cx - hx, cy - hy), grid["size_x"], grid["size_y"],
        linewidth=2.0,
        edgecolor="magenta",
        facecolor="none",
        linestyle="--",
        label="Grid box (x-y proj.)",
        zorder=4,
    )
    ax.add_patch(box_rect)

    # Mark the box centre
    ax.plot(
        cx, cy,
        marker="+", markersize=12,
        color="magenta", markeredgewidth=2,
        label=f"Center ({cx:.1f}, {cy:.1f})",
        zorder=5,
    )

    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.set_title(
        f"Grid box projection — {job_name}\nChain {chain_id}",
        fontsize=12, fontweight="bold",
    )
    ax.legend(fontsize=8, loc="best", framealpha=0.8)
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(alpha=0.3)

    fig.tight_layout()

    if save:
        out_path = _interface_plot_dir() / f"{job_name}_grid_box.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        logger.info("Saved grid box plot: %s", out_path)

    return fig


def plot_switch_distances(
    switch_data: dict,
    save: bool = True,
) -> "plt.Figure | None":
    """Grouped bar chart of Switch I and Switch II distances for GTP vs GDP.

    Parameters
    ----------
    switch_data:
        ``{job_name: {"switch1_mean_dist": float, "switch2_mean_dist": float}}``
        as returned by ``interface.measure_switch_distances()``.
    save:
        If True, save to results/plots/interface/switch_distances.png.

    Returns
    -------
    matplotlib.figure.Figure or None.
    """
    if not _check_mpl("plot_switch_distances"):
        return None

    if not switch_data:
        logger.warning("plot_switch_distances: switch_data is empty.")
        return None

    job_names = list(switch_data.keys())
    sw1_dists = [switch_data[j].get("switch1_mean_dist", float("nan")) for j in job_names]
    sw2_dists = [switch_data[j].get("switch2_mean_dist", float("nan")) for j in job_names]

    # Separate GTP / GDP to assign colours, fall back to sequential palette
    colors: list[str] = []
    for jn in job_names:
        jn_lower = jn.lower()
        if "gtp" in jn_lower:
            colors.append(COLOR_GTP)
        elif "gdp" in jn_lower:
            colors.append(COLOR_GDP)
        else:
            colors.append("#55A868")  # green for other

    n = len(job_names)
    x = np.arange(n)
    bar_width = 0.35

    fig, ax = plt.subplots(figsize=(max(7, n * 1.8), 5))

    bars_sw1 = ax.bar(
        x - bar_width / 2, sw1_dists,
        width=bar_width,
        color=colors,
        alpha=0.85,
        label="Switch I",
        zorder=3,
    )
    bars_sw2 = ax.bar(
        x + bar_width / 2, sw2_dists,
        width=bar_width,
        color=colors,
        alpha=0.55,
        edgecolor="black",
        linewidth=0.7,
        label="Switch II",
        zorder=3,
    )

    # Annotate bars
    for bar in list(bars_sw1) + list(bars_sw2):
        h = bar.get_height()
        if not np.isnan(h):
            ax.annotate(
                f"{h:.1f}",
                xy=(bar.get_x() + bar.get_width() / 2, h),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center", va="bottom",
                fontsize=7,
            )

    # Reference line at 10 Å
    ax.axhline(
        y=10.0,
        color="red",
        linestyle="--",
        linewidth=1.2,
        alpha=0.8,
        label="10 Å threshold (in contact)",
        zorder=2,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(job_names, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("Mean nearest-atom distance (Å)")
    ax.set_title(
        "Switch I & II distances to partner chain\n"
        "(lower = loops toward partner = active conformation)",
        fontsize=12, fontweight="bold",
    )
    ax.grid(axis="y", alpha=0.3, zorder=0)

    # Legend: Switch I / Switch II + GTP / GDP colour patches
    leg_handles = [
        mpatches.Patch(facecolor="white", edgecolor="black", alpha=0.9, label="Switch I (solid)"),
        mpatches.Patch(facecolor="white", edgecolor="black", alpha=0.55, label="Switch II (lighter)"),
        mpatches.Patch(color=COLOR_GTP, label="GTP"),
        mpatches.Patch(color=COLOR_GDP, label="GDP"),
        plt.Line2D([0], [0], color="red", linestyle="--", label="10 Å threshold"),
    ]
    ax.legend(handles=leg_handles, fontsize=8, loc="upper right", framealpha=0.8)

    fig.tight_layout()

    if save:
        out_path = _interface_plot_dir() / "switch_distances.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        logger.info("Saved switch distances plot: %s", out_path)

    return fig
