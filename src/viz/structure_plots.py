"""3D protein structure visualizations using py3Dmol.

Reads CIF files produced by AlphaFold 3 and returns interactive py3Dmol views
suitable for Jupyter notebooks, plus static PNG exports for
results/plots/structures/.

Install dependencies:
    pip install py3Dmol biopython matplotlib
"""

import logging
import os
from pathlib import Path

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Optional py3Dmol import
# ---------------------------------------------------------------------------
try:
    import py3Dmol
    _PY3DMOL_AVAILABLE = True
except ImportError:
    _PY3DMOL_AVAILABLE = False
    logger.warning(
        "py3Dmol is not installed. Install it with:\n"
        "    pip install py3Dmol\n"
        "All visualization functions will return None until it is available."
    )

# ---------------------------------------------------------------------------
# Colour palettes
# ---------------------------------------------------------------------------
CHAIN_COLORS = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2"]

PLDDT_COLORS = {
    "very_high": "#0053D6",   # >= 90  dark blue
    "confident": "#65CBF3",   # 70–90  light blue
    "low":       "#FFDB13",   # 50–70  yellow
    "very_low":  "#FF7D45",   # < 50   orange
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _require_py3dmol():
    """Return True if py3Dmol is available, else log and return False."""
    if not _PY3DMOL_AVAILABLE:
        logger.error(
            "py3Dmol is required for this function. "
            "Install it with: pip install py3Dmol"
        )
        return False
    return True


def _read_cif_string(cif_path: str) -> str:
    """Return the raw CIF file contents as a string."""
    with open(cif_path, "r") as fh:
        return fh.read()


def _get_chains(cif_path: str) -> list:
    """Return a sorted list of unique chain IDs in a CIF file via Biopython."""
    try:
        from Bio.PDB import MMCIFParser
    except ImportError:
        logger.warning("Biopython not installed; chain parsing skipped.")
        return []

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("struct", cif_path)
    chains = []
    seen = set()
    for model in structure:
        for chain in model:
            if chain.id not in seen:
                chains.append(chain.id)
                seen.add(chain.id)
    return sorted(chains)


def _get_residues_with_bfactor(cif_path: str) -> dict:
    """Return {chain_id: [(res_seq, bfactor), ...]} for CA atoms via Biopython."""
    try:
        from Bio.PDB import MMCIFParser
    except ImportError:
        logger.warning("Biopython not installed; residue/bfactor parsing skipped.")
        return {}

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("struct", cif_path)
    result = {}
    for model in structure:
        for chain in model:
            residues = []
            for residue in chain:
                for atom in residue:
                    if atom.name == "CA":
                        residues.append(
                            (residue.id[1], atom.get_bfactor())
                        )
                        break
            if residues:
                result[chain.id] = residues
    return result


def _default_output_dir() -> Path:
    """Resolve the AF3 output directory (AF3_OUTPUT_DIR env var or ../output)."""
    return Path(
        os.environ.get(
            "AF3_OUTPUT_DIR",
            str(Path(__file__).resolve().parent.parent.parent / "output"),
        )
    )


def _structures_plot_dir() -> Path:
    """Return the results/plots/structures/ directory, creating it if needed."""
    repo_root = Path(__file__).resolve().parent.parent.parent
    d = repo_root / "results" / "plots" / "structures"
    d.mkdir(parents=True, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def view_complex_by_chain(cif_path: str, width: int = 800, height: int = 600):
    """Visualize a protein complex coloured by chain.

    Each chain is assigned a distinct colour from the CHAIN_COLORS palette
    (cycles if more than 5 chains). The structure is shown as a cartoon
    representation.

    Parameters
    ----------
    cif_path:
        Path to an AF3 mmCIF file.
    width, height:
        Dimensions of the viewer in pixels.

    Returns
    -------
    py3Dmol.view or None
        Call ``view.show()`` inside a Jupyter notebook to render.

    Examples
    --------
    >>> from src.viz.structure_plots import view_complex_by_chain
    >>> view = view_complex_by_chain("output/my_job/my_job_model.cif")
    >>> view.show()
    """
    if not _require_py3dmol():
        return None

    cif_string = _read_cif_string(cif_path)
    chains = _get_chains(cif_path)

    view = py3Dmol.view(width=width, height=height)
    view.addModel(cif_string, "cif")

    if chains:
        for idx, chain_id in enumerate(chains):
            color = CHAIN_COLORS[idx % len(CHAIN_COLORS)]
            view.setStyle(
                {"chain": chain_id},
                {"cartoon": {"color": color}},
            )
        # Add chain labels at the centroid of each chain's CA atoms
        try:
            residues_by_chain = _get_residues_with_bfactor(cif_path)
            from Bio.PDB import MMCIFParser
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure("struct", cif_path)
            for idx, chain_id in enumerate(chains):
                color = CHAIN_COLORS[idx % len(CHAIN_COLORS)]
                ca_coords = []
                for model in structure:
                    if chain_id in [c.id for c in model]:
                        for residue in model[chain_id]:
                            if "CA" in residue:
                                ca_coords.append(residue["CA"].get_vector().get_array().tolist())
                if ca_coords:
                    cx = sum(c[0] for c in ca_coords) / len(ca_coords)
                    cy = sum(c[1] for c in ca_coords) / len(ca_coords)
                    cz = sum(c[2] for c in ca_coords) / len(ca_coords)
                    view.addLabel(
                        f"Chain {chain_id}",
                        {
                            "position": {"x": cx, "y": cy, "z": cz},
                            "backgroundColor": color,
                            "fontColor": "white",
                            "fontSize": 14,
                            "backgroundOpacity": 0.8,
                        },
                    )
        except Exception as exc:  # noqa: BLE001
            logger.debug("Chain label placement skipped: %s", exc)
    else:
        # Fallback: colour by chain without explicit IDs
        view.setStyle({}, {"cartoon": {"colorscheme": "chain"}})

    view.zoomTo()
    return view


def view_plddt(cif_path: str, width: int = 800, height: int = 600):
    """Visualize a structure coloured by pLDDT confidence score.

    pLDDT is stored in the B-factor column of AF3 CIF files.  The standard
    AF3 colour scheme is used:

    * >= 90  — dark blue  (#0053D6, very high confidence)
    * 70–90  — light blue (#65CBF3, confident)
    * 50–70  — yellow     (#FFDB13, low confidence)
    * <  50  — orange     (#FF7D45, very low confidence)

    Parameters
    ----------
    cif_path:
        Path to an AF3 mmCIF file.
    width, height:
        Viewer dimensions in pixels.

    Returns
    -------
    py3Dmol.view or None

    Examples
    --------
    >>> from src.viz.structure_plots import view_plddt
    >>> view = view_plddt("output/my_job/my_job_model.cif")
    >>> view.show()
    """
    if not _require_py3dmol():
        return None

    cif_string = _read_cif_string(cif_path)

    # Gather per-residue pLDDT (from B-factor of CA atoms)
    residues_by_chain = _get_residues_with_bfactor(cif_path)

    view = py3Dmol.view(width=width, height=height)
    view.addModel(cif_string, "cif")

    if residues_by_chain:
        for chain_id, res_list in residues_by_chain.items():
            for res_seq, bfactor in res_list:
                plddt = bfactor
                if plddt >= 90:
                    color = PLDDT_COLORS["very_high"]
                elif plddt >= 70:
                    color = PLDDT_COLORS["confident"]
                elif plddt >= 50:
                    color = PLDDT_COLORS["low"]
                else:
                    color = PLDDT_COLORS["very_low"]

                view.setStyle(
                    {"chain": chain_id, "resi": res_seq},
                    {"cartoon": {"color": color}},
                )
    else:
        # If Biopython is unavailable fall back to built-in b-factor coloring
        view.setStyle(
            {},
            {"cartoon": {"colorscheme": {"prop": "b", "gradient": "roygb", "min": 50, "max": 90}}},
        )

    view.zoomTo()
    return view


def view_binding_site(
    cif_path: str,
    chain_a_id: str,
    chain_b_id: str,
    contact_residues_a: list,
    contact_residues_b: list,
    width: int = 800,
    height: int = 600,
):
    """Visualize the interface between two chains.

    The full complex is shown as a faded cartoon (opacity 0.3).  Contact
    residues on each chain are highlighted as opaque sticks in contrasting
    colours, and the camera zooms to the interface region.

    Parameters
    ----------
    cif_path:
        Path to an AF3 mmCIF file.
    chain_a_id:
        Chain ID of the first partner (e.g. ``"A"``).
    chain_b_id:
        Chain ID of the second partner (e.g. ``"B"``).
    contact_residues_a:
        List of residue sequence numbers on chain A that form contacts.
    contact_residues_b:
        List of residue sequence numbers on chain B that form contacts.
    width, height:
        Viewer dimensions in pixels.

    Returns
    -------
    py3Dmol.view or None

    Examples
    --------
    >>> view = view_binding_site(
    ...     "output/rab5_effector/model.cif",
    ...     chain_a_id="A", chain_b_id="B",
    ...     contact_residues_a=[40, 41, 42, 70, 71],
    ...     contact_residues_b=[10, 11, 12],
    ... )
    >>> view.show()
    """
    if not _require_py3dmol():
        return None

    cif_string = _read_cif_string(cif_path)

    view = py3Dmol.view(width=width, height=height)
    view.addModel(cif_string, "cif")

    # Whole complex as faded cartoon
    view.setStyle({}, {"cartoon": {"color": "lightgrey", "opacity": 0.3}})

    # Contact residues as opaque sticks
    color_a = "#4C72B0"  # blue
    color_b = "#DD8452"  # orange

    if contact_residues_a:
        resi_str_a = ",".join(str(r) for r in contact_residues_a)
        view.setStyle(
            {"chain": chain_a_id, "resi": resi_str_a},
            {"stick": {"color": color_a, "opacity": 1.0}},
        )
    if contact_residues_b:
        resi_str_b = ",".join(str(r) for r in contact_residues_b)
        view.setStyle(
            {"chain": chain_b_id, "resi": resi_str_b},
            {"stick": {"color": color_b, "opacity": 1.0}},
        )

    # Zoom to interface: build a selection that covers both contact sets
    contact_sel = []
    if contact_residues_a:
        contact_sel.append(
            {"chain": chain_a_id, "resi": contact_residues_a}
        )
    if contact_residues_b:
        contact_sel.append(
            {"chain": chain_b_id, "resi": contact_residues_b}
        )

    if contact_sel:
        # py3Dmol zoomTo accepts a selection dict; use the first chain's
        # contact residues as anchor and let the viewer frame both.
        view.zoomTo({"chain": chain_a_id, "resi": contact_residues_a} if contact_residues_a else {})
    else:
        view.zoomTo()

    return view


def view_gtp_gdp_overlay(
    gtp_cif_path: str,
    gdp_cif_path: str,
    rab_chain_id: str = "A",
    switch1_range: tuple = (40, 52),
    switch2_range: tuple = (69, 82),
    width: int = 800,
    height: int = 600,
):
    """Overlay GTP-bound and GDP-bound structures highlighting switch regions.

    Both structures are added to the same py3Dmol viewer in their original
    coordinate frames.  The GTP structure is shown in blue, the GDP structure
    in orange.  Switch I and Switch II residues are highlighted as thick sticks.

    .. note::
        py3Dmol does **not** perform structural alignment.  If the two
        structures are in different coordinate frames (e.g. from separate AF3
        runs without common reference), pre-align them with Biopython or
        MDAnalysis before calling this function.  A simple approach:

        .. code-block:: python

            from Bio.PDB import Superimposer, MMCIFParser, MMCIFIO
            # Align gdp_structure onto gtp_structure using backbone atoms,
            # then save the transformed GDP structure to a new CIF file.

    Parameters
    ----------
    gtp_cif_path:
        Path to the GTP-bound AF3 CIF file.
    gdp_cif_path:
        Path to the GDP-bound AF3 CIF file.
    rab_chain_id:
        Chain ID of the Rab GTPase (default ``"A"``).
    switch1_range:
        (start, end) residue numbers for Switch I (inclusive).
    switch2_range:
        (start, end) residue numbers for Switch II (inclusive).
    width, height:
        Viewer dimensions in pixels.

    Returns
    -------
    py3Dmol.view or None

    Examples
    --------
    >>> view = view_gtp_gdp_overlay(
    ...     "output/rab5_gtp/model.cif",
    ...     "output/rab5_gdp/model.cif",
    ...     rab_chain_id="A",
    ...     switch1_range=(40, 52),
    ...     switch2_range=(69, 82),
    ... )
    >>> view.show()
    """
    if not _require_py3dmol():
        return None

    gtp_string = _read_cif_string(gtp_cif_path)
    gdp_string = _read_cif_string(gdp_cif_path)

    sw1_start, sw1_end = switch1_range
    sw2_start, sw2_end = switch2_range
    sw1_resi = list(range(sw1_start, sw1_end + 1))
    sw2_resi = list(range(sw2_start, sw2_end + 1))
    sw_resi = sw1_resi + sw2_resi
    sw_resi_str = ",".join(str(r) for r in sw_resi)

    color_gtp = "#4C72B0"   # blue
    color_gdp = "#DD8452"   # orange
    color_sw1 = "#55A868"   # green — Switch I
    color_sw2 = "#C44E52"   # red   — Switch II

    view = py3Dmol.view(width=width, height=height)

    # ------ GTP model (model index 0) ------
    view.addModel(gtp_string, "cif")
    # Whole chain as thin cartoon
    view.setStyle(
        {"model": 0, "chain": rab_chain_id},
        {"cartoon": {"color": color_gtp, "thickness": 0.4}},
    )
    # Switch I thick sticks
    view.setStyle(
        {"model": 0, "chain": rab_chain_id, "resi": ",".join(str(r) for r in sw1_resi)},
        {"stick": {"color": color_sw1, "radius": 0.35}},
    )
    # Switch II thick sticks
    view.setStyle(
        {"model": 0, "chain": rab_chain_id, "resi": ",".join(str(r) for r in sw2_resi)},
        {"stick": {"color": color_sw2, "radius": 0.35}},
    )

    # ------ GDP model (model index 1) ------
    view.addModel(gdp_string, "cif")
    view.setStyle(
        {"model": 1, "chain": rab_chain_id},
        {"cartoon": {"color": color_gdp, "thickness": 0.4}},
    )
    view.setStyle(
        {"model": 1, "chain": rab_chain_id, "resi": ",".join(str(r) for r in sw1_resi)},
        {"stick": {"color": color_sw1, "radius": 0.35}},
    )
    view.setStyle(
        {"model": 1, "chain": rab_chain_id, "resi": ",".join(str(r) for r in sw2_resi)},
        {"stick": {"color": color_sw2, "radius": 0.35}},
    )

    # Labels for switch regions (based on GTP model)
    view.addLabel(
        "Switch I",
        {
            "position": {"resi": sw1_start, "chain": rab_chain_id, "model": 0},
            "backgroundColor": color_sw1,
            "fontColor": "white",
            "fontSize": 12,
        },
    )
    view.addLabel(
        "Switch II",
        {
            "position": {"resi": sw2_start, "chain": rab_chain_id, "model": 0},
            "backgroundColor": color_sw2,
            "fontColor": "white",
            "fontSize": 12,
        },
    )

    view.zoomTo({"chain": rab_chain_id, "resi": sw_resi_str})
    return view


def view_docking_grid(
    receptor_cif_path: str,
    chain_id: str,
    grid: dict,
    width: int = 800,
    height: int = 600,
):
    """Visualize receptor surface with docking grid box overlay.

    The receptor is shown as a semi-transparent surface with a cartoon
    underneath.  The docking search box is drawn as a wireframe using
    py3Dmol's ``addBox()`` if available, otherwise as 12 cylinder edges.

    Parameters
    ----------
    receptor_cif_path:
        Path to the receptor AF3 CIF file.
    chain_id:
        Chain ID of the receptor to display.
    grid:
        Dictionary with keys:
        ``center_x``, ``center_y``, ``center_z``,
        ``size_x``, ``size_y``, ``size_z`` (all in Angstroms).
    width, height:
        Viewer dimensions in pixels.

    Returns
    -------
    py3Dmol.view or None

    Examples
    --------
    >>> grid = {
    ...     "center_x": 10.0, "center_y": 20.0, "center_z": 5.0,
    ...     "size_x": 20.0, "size_y": 20.0, "size_z": 20.0,
    ... }
    >>> view = view_docking_grid("output/my_receptor/model.cif", "A", grid)
    >>> view.show()
    """
    if not _require_py3dmol():
        return None

    cif_string = _read_cif_string(receptor_cif_path)

    view = py3Dmol.view(width=width, height=height)
    view.addModel(cif_string, "cif")

    # Receptor: semi-transparent surface + cartoon
    chain_sel = {"chain": chain_id}
    view.addSurface(
        py3Dmol.SAS,
        {"opacity": 0.5, "color": "lightblue"},
        chain_sel,
    )
    view.setStyle(chain_sel, {"cartoon": {"color": "lightblue", "opacity": 0.8}})

    # Draw docking box
    cx = grid["center_x"]
    cy = grid["center_y"]
    cz = grid["center_z"]
    hx = grid["size_x"] / 2.0
    hy = grid["size_y"] / 2.0
    hz = grid["size_z"] / 2.0

    box_spec = {
        "center": {"x": cx, "y": cy, "z": cz},
        "dimensions": {"w": grid["size_x"], "h": grid["size_y"], "d": grid["size_z"]},
        "color": "magenta",
        "opacity": 0.9,
        "wireframe": True,
    }

    # Try addBox first; fall back to 12 cylinder edges
    try:
        view.addBox(box_spec)
    except AttributeError:
        logger.debug("py3Dmol.addBox not available; drawing 12 cylinder edges instead.")
        _draw_box_edges(view, cx, cy, cz, hx, hy, hz)

    view.zoomTo()
    return view


def _draw_box_edges(view, cx, cy, cz, hx, hy, hz):
    """Draw 12 edges of an axis-aligned box using thin cylinders."""
    corners = [
        (cx - hx, cy - hy, cz - hz),
        (cx + hx, cy - hy, cz - hz),
        (cx + hx, cy + hy, cz - hz),
        (cx - hx, cy + hy, cz - hz),
        (cx - hx, cy - hy, cz + hz),
        (cx + hx, cy - hy, cz + hz),
        (cx + hx, cy + hy, cz + hz),
        (cx - hx, cy + hy, cz + hz),
    ]
    edges = [
        (0, 1), (1, 2), (2, 3), (3, 0),  # bottom face
        (4, 5), (5, 6), (6, 7), (7, 4),  # top face
        (0, 4), (1, 5), (2, 6), (3, 7),  # vertical pillars
    ]
    for i, j in edges:
        x1, y1, z1 = corners[i]
        x2, y2, z2 = corners[j]
        view.addCylinder(
            {
                "start": {"x": x1, "y": y1, "z": z1},
                "end":   {"x": x2, "y": y2, "z": z2},
                "radius": 0.15,
                "color": "magenta",
                "fromCap": 1,
                "toCap": 1,
            }
        )


def view_docked_pose(
    receptor_pdbqt_path: str,
    ligand_pose_path: str,
    grid: dict,
    width: int = 800,
    height: int = 600,
):
    """Visualize a docked receptor–ligand pose with grid box.

    The receptor is shown as a semi-transparent light-gray surface.
    The ligand is shown as CPK-colored sticks.  The docking grid box is
    drawn as a wireframe.

    Parameters
    ----------
    receptor_pdbqt_path:
        Path to the receptor file (PDBQT or PDB).
    ligand_pose_path:
        Path to the docked ligand pose file (PDBQT, SDF, MOL2, or PDB).
    grid:
        Dictionary with keys:
        ``center_x``, ``center_y``, ``center_z``,
        ``size_x``, ``size_y``, ``size_z``.
    width, height:
        Viewer dimensions in pixels.

    Returns
    -------
    py3Dmol.view or None

    Examples
    --------
    >>> grid = {
    ...     "center_x": 10.0, "center_y": 20.0, "center_z": 5.0,
    ...     "size_x": 20.0, "size_y": 20.0, "size_z": 20.0,
    ... }
    >>> view = view_docked_pose(
    ...     "docking/receptor.pdbqt",
    ...     "docking/best_pose.pdbqt",
    ...     grid,
    ... )
    >>> view.show()
    """
    if not _require_py3dmol():
        return None

    # Determine format from extension
    rec_ext = Path(receptor_pdbqt_path).suffix.lower().lstrip(".")
    lig_ext = Path(ligand_pose_path).suffix.lower().lstrip(".")
    # py3Dmol accepts "pdb" for both pdb and pdbqt
    rec_fmt = "pdb"
    lig_fmt = lig_ext if lig_ext in ("sdf", "mol2", "xyz") else "pdb"

    rec_string = _read_cif_string(receptor_pdbqt_path)
    lig_string = _read_cif_string(ligand_pose_path)

    view = py3Dmol.view(width=width, height=height)

    # Receptor: semi-transparent surface + light gray cartoon
    view.addModel(rec_string, rec_fmt)
    view.addSurface(
        py3Dmol.SAS,
        {"opacity": 0.4, "color": "lightgrey"},
        {"model": 0},
    )
    view.setStyle({"model": 0}, {"cartoon": {"color": "lightgrey", "opacity": 0.6}})

    # Ligand: CPK-colored sticks
    view.addModel(lig_string, lig_fmt)
    view.setStyle({"model": 1}, {"stick": {"colorscheme": "default"}})

    # Grid box
    cx = grid["center_x"]
    cy = grid["center_y"]
    cz = grid["center_z"]
    hx = grid["size_x"] / 2.0
    hy = grid["size_y"] / 2.0
    hz = grid["size_z"] / 2.0

    box_spec = {
        "center": {"x": cx, "y": cy, "z": cz},
        "dimensions": {"w": grid["size_x"], "h": grid["size_y"], "d": grid["size_z"]},
        "color": "magenta",
        "opacity": 0.9,
        "wireframe": True,
    }
    try:
        view.addBox(box_spec)
    except AttributeError:
        _draw_box_edges(view, cx, cy, cz, hx, hy, hz)

    view.zoomTo({"model": 1})
    return view


def save_view_png(view, output_path: str) -> None:
    """Save a py3Dmol view to a PNG file.

    Attempts to use ``view.png()`` (available in recent py3Dmol builds running
    in a headless environment).  If that is not available — which is common
    in plain Jupyter / IPython contexts — a placeholder PNG is written via
    matplotlib so that the ``results/plots/structures/`` directory is always
    populated.

    Parameters
    ----------
    view:
        A ``py3Dmol.view`` object.
    output_path:
        Destination file path (should end in ``.png``).

    Examples
    --------
    >>> from src.viz.structure_plots import view_plddt, save_view_png
    >>> v = view_plddt("output/my_job/model.cif")
    >>> save_view_png(v, "results/plots/structures/my_job/plddt.png")
    """
    output_path = str(output_path)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    if view is None:
        _write_placeholder_png(output_path, "No py3Dmol view available")
        return

    # Try view.png() — may work in some environments
    try:
        png_data = view.png()
        if png_data:
            import base64
            # png() may return a data-URI string: "data:image/png;base64,..."
            if isinstance(png_data, str) and png_data.startswith("data:"):
                b64 = png_data.split(",", 1)[1]
                raw = base64.b64decode(b64)
            elif isinstance(png_data, (bytes, bytearray)):
                raw = bytes(png_data)
            else:
                raw = base64.b64decode(png_data)

            with open(output_path, "wb") as fh:
                fh.write(raw)
            logger.info("Saved PNG via view.png(): %s", output_path)
            return
    except Exception as exc:  # noqa: BLE001
        logger.warning(
            "view.png() failed (%s); writing placeholder PNG instead.", exc
        )

    _write_placeholder_png(
        output_path,
        "3D view: open notebook for\ninteractive visualization",
    )


def _write_placeholder_png(output_path: str, message: str) -> None:
    """Write a simple text-on-white matplotlib figure as a placeholder PNG."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning(
            "matplotlib is not installed; cannot write placeholder PNG for %s. "
            "Install it with: pip install matplotlib",
            output_path,
        )
        return

    fig, ax = plt.subplots(figsize=(6, 3))
    ax.text(
        0.5, 0.5, message,
        ha="center", va="center",
        fontsize=14, wrap=True,
        transform=ax.transAxes,
    )
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(output_path, dpi=96, bbox_inches="tight")
    plt.close(fig)
    logger.warning(
        "Saved placeholder PNG (py3Dmol PNG export not available): %s",
        output_path,
    )


def make_structure_report(job_name: str, output_dir: str = None) -> None:
    """Generate standard structure visualizations for an AF3 job.

    Finds the top-level CIF model file for *job_name*, creates complex-by-chain
    and pLDDT views, and saves PNG exports to
    ``results/plots/structures/<job_name>/``.

    Views that require additional inputs (contact residues, docking grid) are
    skipped automatically.

    Parameters
    ----------
    job_name:
        Name of the AF3 job directory (e.g. ``"rab5_gtp_effector"``).
    output_dir:
        Root AF3 output directory.  Resolved in this order:
        1. *output_dir* argument (if not None)
        2. ``AF3_OUTPUT_DIR`` environment variable
        3. ``../output`` relative to the repo root

    Examples
    --------
    >>> from src.viz.structure_plots import make_structure_report
    >>> make_structure_report("rab5_gtp_effector")
    # Creates:
    #   results/plots/structures/rab5_gtp_effector/complex_by_chain.png
    #   results/plots/structures/rab5_gtp_effector/plddt.png
    """
    if output_dir is None:
        base = _default_output_dir()
    else:
        base = Path(output_dir)

    job_dir = base / job_name

    # Locate the top-level CIF model (AF3 naming: <job_name>_model.cif)
    cif_candidates = list(job_dir.glob(f"{job_name}_model.cif"))
    if not cif_candidates:
        # Fallback: any CIF directly in the job directory
        cif_candidates = list(job_dir.glob("*.cif"))
        # Prefer files not inside seed-*_sample-* subdirs
        cif_candidates = [c for c in cif_candidates if c.parent == job_dir]

    if not cif_candidates:
        logger.error(
            "No CIF file found for job '%s' in %s", job_name, job_dir
        )
        return

    cif_path = str(cif_candidates[0])
    logger.info("make_structure_report: using CIF %s", cif_path)

    out_dir = _structures_plot_dir() / job_name
    out_dir.mkdir(parents=True, exist_ok=True)

    views = {
        "complex_by_chain": view_complex_by_chain,
        "plddt":             view_plddt,
    }

    for name, fn in views.items():
        try:
            v = fn(cif_path)
            save_view_png(v, str(out_dir / f"{name}.png"))
        except Exception as exc:  # noqa: BLE001
            logger.error("Failed to generate %s view for %s: %s", name, job_name, exc)

    logger.info(
        "make_structure_report complete. PNGs saved to %s", out_dir
    )
