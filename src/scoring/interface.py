"""Compute interface metrics from an AF3 CIF structure file using Biopython."""
import logging
import warnings
from pathlib import Path
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Switch-region definitions for RAB GTPases (1-indexed sequence positions)
# ---------------------------------------------------------------------------
SWITCH1_RANGE = (40, 52)   # inclusive
SWITCH2_RANGE = (69, 82)   # inclusive


def _load_structure(cif_path: str):
    """Parse a CIF file and return the first model."""
    from Bio.PDB import MMCIFParser  # type: ignore

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("struct", cif_path)
    return next(iter(structure))  # first model


def _get_chain(model, chain_id: str):
    """Return chain object or None, with a warning if missing."""
    try:
        return model[chain_id]
    except KeyError:
        logger.warning("Chain '%s' not found in structure.", chain_id)
        return None


def _chain_atoms(chain) -> list:
    """Return a flat list of all Atom objects in a chain."""
    return list(chain.get_atoms())


def _residue_index(residue) -> int:
    """Return 0-based sequential index from residue sequence number."""
    # Biopython residue id is (hetflag, seq_id, icode); seq_id is 1-based
    return residue.get_id()[1] - 1


# ---------------------------------------------------------------------------
# BSA estimation
# ---------------------------------------------------------------------------

def _sasa_sphere_estimate(atoms, probe_radius: float = 1.4) -> float:
    """
    Rough SASA estimate via Lee-Richards-style shrink-wrap approximation.
    Each atom contributes 4π(r+probe)² reduced by neighbor overlaps.
    This is only used when DSSP is unavailable.
    """
    # Van-der-Waals radii by element (Å)
    VDW = {"C": 1.7, "N": 1.55, "O": 1.52, "S": 1.80,
           "P": 1.80, "H": 1.20}
    DEFAULT_VDW = 1.7

    coords = np.array([a.get_vector().get_array() for a in atoms])
    radii = np.array([VDW.get(a.element.strip().upper() if a.element else "C",
                              DEFAULT_VDW) + probe_radius
                      for a in atoms])

    total = 0.0
    n = len(coords)
    if n == 0:
        return 0.0

    for i in range(n):
        sphere_area = 4.0 * np.pi * radii[i] ** 2
        # fraction buried by neighbouring atoms
        diffs = coords - coords[i]
        dists = np.linalg.norm(diffs, axis=1)
        dists[i] = np.inf  # exclude self
        overlapping = dists < (radii + radii[i])
        # Each overlap reduces exposed area proportionally
        buried_frac = 0.0
        for j in np.where(overlapping)[0]:
            d = dists[j]
            ri, rj = radii[i], radii[j]
            # fractional cap buried by atom j
            h = (ri + rj - d) / (2 * ri)
            buried_frac += max(0.0, min(h, 1.0))
        exposed = sphere_area * max(0.0, 1.0 - buried_frac)
        total += exposed
    return total


def _compute_bsa(chain_a, chain_b) -> float:
    """
    Approximate BSA (Å²) = (SASA_A_free + SASA_B_free) − SASA_complex.
    Uses Biopython DSSP when available, else a sphere-overlap estimate.
    """
    atoms_a = _chain_atoms(chain_a)
    atoms_b = _chain_atoms(chain_b)
    atoms_ab = atoms_a + atoms_b

    try:
        from Bio.PDB.SASA import ShrakeRupley  # type: ignore

        sr = ShrakeRupley()

        # Build minimal fake structures for SR computation
        from Bio.PDB import Structure, Model, Chain  # type: ignore

        def _mini_struct(chain_obj, struct_id: str):
            s = Structure.Structure(struct_id)
            m = Model.Model(0)
            s.add(m)
            import copy
            c = copy.deepcopy(chain_obj)
            m.add(c)
            return s

        s_a = _mini_struct(chain_a, "A")
        s_b = _mini_struct(chain_b, "B")

        import copy
        from Bio.PDB import Structure, Model  # type: ignore

        s_ab = Structure.Structure("AB")
        m_ab = Model.Model(0)
        s_ab.add(m_ab)
        m_ab.add(copy.deepcopy(chain_a))
        m_ab.add(copy.deepcopy(chain_b))

        sr.compute(s_a, level="A")
        sr.compute(s_b, level="A")
        sr.compute(s_ab, level="A")

        def _total_sasa(struct):
            return sum(a.sasa for a in struct.get_atoms())

        sasa_a = _total_sasa(s_a)
        sasa_b = _total_sasa(s_b)
        sasa_ab = _total_sasa(s_ab)
        return max(0.0, (sasa_a + sasa_b) - sasa_ab)

    except Exception as exc:
        logger.debug("ShrakeRupley unavailable (%s); using sphere-overlap estimate.", exc)
        sasa_a = _sasa_sphere_estimate(atoms_a)
        sasa_b = _sasa_sphere_estimate(atoms_b)
        sasa_ab = _sasa_sphere_estimate(atoms_ab)
        return max(0.0, (sasa_a + sasa_b) - sasa_ab)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute_interface(
    cif_path: str,
    chain_a_id: str,
    chain_b_id: str,
    cutoff: float = 5.0,
) -> dict:
    """Compute interface metrics between two chains in an AF3 CIF file.

    Parameters
    ----------
    cif_path:
        Path to the AF3 CIF model file.
    chain_a_id:
        Chain identifier for chain A (e.g. "A").
    chain_b_id:
        Chain identifier for chain B (e.g. "B").
    cutoff:
        Distance cutoff in Å for contact detection (default 5.0 Å).

    Returns
    -------
    dict with keys:
        bsa (float): buried surface area in Å²
        contact_count (int): number of residue pairs within cutoff
        contact_residues_a (list[int]): 0-indexed residue positions in chain A
        contact_residues_b (list[int]): 0-indexed residue positions in chain B
    """
    from Bio.PDB import NeighborSearch  # type: ignore

    model = _load_structure(cif_path)
    chain_a = _get_chain(model, chain_a_id)
    chain_b = _get_chain(model, chain_b_id)

    if chain_a is None or chain_b is None:
        logger.warning(
            "compute_interface: missing chain(s) %s/%s — returning empty dict.",
            chain_a_id, chain_b_id,
        )
        return {}

    # NeighborSearch over all atoms in chain B; query with atoms from chain A
    atoms_b = _chain_atoms(chain_b)
    ns = NeighborSearch(atoms_b)

    contact_pairs: set[tuple[int, int]] = set()

    for atom_a in _chain_atoms(chain_a):
        res_a = atom_a.get_parent()
        idx_a = _residue_index(res_a)
        nearby = ns.search(atom_a.get_vector().get_array(), cutoff, level="A")
        for atom_b in nearby:
            res_b = atom_b.get_parent()
            idx_b = _residue_index(res_b)
            contact_pairs.add((idx_a, idx_b))

    contact_residues_a = sorted({p[0] for p in contact_pairs})
    contact_residues_b = sorted({p[1] for p in contact_pairs})

    bsa = _compute_bsa(chain_a, chain_b)

    return {
        "bsa": bsa,
        "contact_count": len(contact_pairs),
        "contact_residues_a": contact_residues_a,
        "contact_residues_b": contact_residues_b,
    }


def get_contact_residues(
    cif_path: str,
    chain_a_id: str,
    chain_b_id: str,
    cutoff: float = 5.0,
) -> tuple[list[int], list[int]]:
    """Convenience wrapper returning (residues_in_a, residues_in_b).

    Parameters
    ----------
    cif_path:
        Path to the AF3 CIF model file.
    chain_a_id:
        Chain identifier for chain A.
    chain_b_id:
        Chain identifier for chain B.
    cutoff:
        Distance cutoff in Å (default 5.0 Å).

    Returns
    -------
    Tuple of (contact_residues_a, contact_residues_b) as 0-indexed lists.
    """
    result = compute_interface(cif_path, chain_a_id, chain_b_id, cutoff)
    if not result:
        return [], []
    return result["contact_residues_a"], result["contact_residues_b"]


def measure_switch_distances(
    cif_path: str,
    rab_chain_id: str,
    partner_chain_id: str,
) -> dict:
    """Measure mean distance from RAB GTPase switch regions to partner chain.

    Switch I  ~ residues 40–52 (1-indexed), Switch II ~ residues 69–82.
    Lower mean distance indicates loops pointing toward partner
    (active / GTP-like conformation).

    Parameters
    ----------
    cif_path:
        Path to the AF3 CIF model file.
    rab_chain_id:
        Chain ID for the RAB GTPase chain.
    partner_chain_id:
        Chain ID for the partner (effector) chain.

    Returns
    -------
    dict with keys:
        switch1_mean_dist (float): mean nearest-atom distance for Switch I (Å)
        switch2_mean_dist (float): mean nearest-atom distance for Switch II (Å)
    """
    model = _load_structure(cif_path)
    rab_chain = _get_chain(model, rab_chain_id)
    partner_chain = _get_chain(model, partner_chain_id)

    if rab_chain is None or partner_chain is None:
        logger.warning(
            "measure_switch_distances: missing chain(s) %s/%s.",
            rab_chain_id, partner_chain_id,
        )
        return {"switch1_mean_dist": float("nan"), "switch2_mean_dist": float("nan")}

    partner_atoms = _chain_atoms(partner_chain)
    if not partner_atoms:
        logger.warning("Partner chain %s has no atoms.", partner_chain_id)
        return {"switch1_mean_dist": float("nan"), "switch2_mean_dist": float("nan")}

    partner_coords = np.array([a.get_vector().get_array() for a in partner_atoms])

    def _mean_nearest_dist(seq_start: int, seq_end: int) -> float:
        """Mean over switch residues of the minimum distance to any partner atom."""
        per_residue: list[float] = []
        for residue in rab_chain.get_residues():
            seq_id = residue.get_id()[1]  # 1-indexed
            if seq_start <= seq_id <= seq_end:
                res_atoms = list(residue.get_atoms())
                if not res_atoms:
                    continue
                atom_dists = []
                for atom in res_atoms:
                    coord = atom.get_vector().get_array()
                    dists = np.linalg.norm(partner_coords - coord, axis=1)
                    atom_dists.append(float(dists.min()))
                per_residue.append(min(atom_dists))  # closest atom in residue
        return float(np.mean(per_residue)) if per_residue else float("nan")

    switch1_mean = _mean_nearest_dist(*SWITCH1_RANGE)
    switch2_mean = _mean_nearest_dist(*SWITCH2_RANGE)

    return {
        "switch1_mean_dist": switch1_mean,
        "switch2_mean_dist": switch2_mean,
    }


if __name__ == "__main__":
    import sys
    import os

    logging.basicConfig(level=logging.INFO, stream=sys.stderr)

    # Quick self-test: look for any CIF in the output directory
    output_dir = os.environ.get(
        "AF3_OUTPUT_DIR",
        str(Path(__file__).resolve().parent.parent.parent / "output"),
    )
    cif_files = list(Path(output_dir).glob("**/*.cif"))
    if not cif_files:
        print("No CIF files found for self-test — skipping.")
        sys.exit(0)

    test_cif = str(cif_files[0])
    print(f"Testing with: {test_cif}")

    from Bio.PDB import MMCIFParser  # type: ignore
    parser = MMCIFParser(QUIET=True)
    struct = parser.get_structure("s", test_cif)
    model = next(iter(struct))
    chain_ids = [c.get_id() for c in model.get_chains()]
    print(f"Chains in file: {chain_ids}")

    if len(chain_ids) >= 2:
        result = compute_interface(test_cif, chain_ids[0], chain_ids[1])
        print("compute_interface result:", result)
        a_res, b_res = get_contact_residues(test_cif, chain_ids[0], chain_ids[1])
        print(f"Contact residues A: {a_res[:10]}{'...' if len(a_res) > 10 else ''}")
        print(f"Contact residues B: {b_res[:10]}{'...' if len(b_res) > 10 else ''}")
        sw = measure_switch_distances(test_cif, chain_ids[0], chain_ids[1])
        print("Switch distances:", sw)
    else:
        print("Only one chain found; skipping multi-chain tests.")
