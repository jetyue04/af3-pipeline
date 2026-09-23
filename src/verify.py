"""
pocket_verify.py

Reusable binding-site / interface verification for AlphaFold3 cofold outputs.
Automates the checks we did by hand in ChimeraX for Rab5A-Rabaptin5 (P-loop
occupancy, switch-region overlap, pairwise ipTM, buried interface area,
GTP-vs-GDP interface comparison) so they can be run across every AF3 job in
the FAM129B / RAB5 / RAB7 / CapZ screen, including the Tier-1 positive
controls (EEA1, RILP) from guide.md.

Expected AF3 output layout per job (matches what you've been generating):

    <job_name>/
        <job_name>_data.json                 # input spec: sequences in chain order
        <job_name>_summary_confidences.json  # top-line ipTM/ptm etc (per ranked-1 sample by default*)
        <job_name>_confidences.json          # full per-sample confidences (pae, chain_pair_iptm, etc)
        ranking_scores.csv                   # seed,sample,ranking_score
        seed-<n>_sample-<m>/
            model.cif

    * NOTE: AF3 write one summary/confidences JSON per sample directory in
      newer releases (seed-*_sample-*/*_summary_confidences.json). This module
      handles both layouts -- see _find_confidence_files().

Dependencies:
    pip install MDAnalysis freesasa pandas numpy --break-system-packages

Typical usage
-------------
    from pocket_verify import score_job, batch_score, compare_nucleotide_states

    # One job
    result = score_job("rab5a_gtp_mg_rabaptin5_complex", chain_a="A", chain_b="B")
    print(result.verdict, result.pairwise_iptm, result.buried_area)

    # Whole screen directory (many job folders as subdirectories)
    df = batch_score("/path/to/af3_outputs", chain_a="A", chain_b="B")
    df.to_csv("results/scores_summary.csv", index=False)

    # GTP vs GDP interface diff (nucleotide-state selectivity check)
    diff = compare_nucleotide_states(
        gtp_job_dir="rab5a_gtp_mg_rabaptin5_complex",
        gdp_job_dir="rab5a_gdp_mg_rabaptin5_complex",
        chain_a="A", chain_b="B",
    )
    print(diff)
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import distances as mda_distances
except ImportError as e:
    raise ImportError(
        "MDAnalysis is required: pip install MDAnalysis --break-system-packages"
    ) from e

try:
    import freesasa
except ImportError:
    freesasa = None  # buried-area calc will raise a clear error if actually called

warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")


# ---------------------------------------------------------------------------
# Thresholds -- guide.md step 7 "verdict table", made explicit and centralized
# so they're tuned in one place instead of copy-pasted across scripts.
# ---------------------------------------------------------------------------

IPTM_STRONG = 0.80
IPTM_WEAK = 0.60
MIN_BURIED_AREA_A2 = 400.0     # below this, treat as likely crystallographic/AF noise, not real
CONSISTENCY_IPTM_STD_MAX = 0.10  # max stdev of ranked-model ipTM to call it "consistent"


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class JobScore:
    job_name: str
    chain_a: str
    chain_b: str
    top_sample_dir: Optional[str] = None
    ranking_score: Optional[float] = None
    global_iptm: Optional[float] = None
    global_ptm: Optional[float] = None
    pairwise_iptm: Optional[float] = None
    has_clash: Optional[bool] = None
    fraction_disordered: Optional[float] = None
    n_interface_residues_a: Optional[int] = None
    n_interface_residues_b: Optional[int] = None
    buried_area_a2: Optional[float] = None
    interface_residues_a: list = field(default_factory=list)
    interface_residues_b: list = field(default_factory=list)
    ranked_iptm_values: list = field(default_factory=list)  # across all samples, for consistency check
    verdict: str = "unscored"
    notes: str = ""

    def as_row(self) -> dict:
        """Flatten to a dict suitable for a pandas DataFrame row."""
        d = self.__dict__.copy()
        d["interface_residues_a"] = ",".join(map(str, self.interface_residues_a))
        d["interface_residues_b"] = ",".join(map(str, self.interface_residues_b))
        d["ranked_iptm_values"] = ",".join(f"{v:.3f}" for v in self.ranked_iptm_values)
        return d


# ---------------------------------------------------------------------------
# File discovery -- handles both AF3 output layouts (job-level vs per-sample
# confidence JSONs), since this has shifted across AF3 releases.
# ---------------------------------------------------------------------------

def _find_confidence_files(job_dir: Path, job_name: str) -> dict:
    """Locate summary_confidences.json, confidences.json, ranking_scores.csv,
    and each seed-*_sample-*/model.cif, regardless of which AF3 layout was used.
    """
    files = {
        "summary": None,
        "full": None,
        "ranking_csv": None,
        "samples": {},  # sample_dir_name -> model.cif path
    }

    # job-level files (older / default layout)
    cand_summary = job_dir / f"{job_name}_summary_confidences.json"
    cand_full = job_dir / f"{job_name}_confidences.json"
    cand_ranking = job_dir / "ranking_scores.csv"
    if cand_summary.exists():
        files["summary"] = cand_summary
    if cand_full.exists():
        files["full"] = cand_full
    if cand_ranking.exists():
        files["ranking_csv"] = cand_ranking

    for sample_dir in sorted(job_dir.glob("seed-*_sample-*")):
        cif = sample_dir / "model.cif"
        if not cif.exists():
            # some layouts name it <job_name>_model.cif inside the sample dir
            alt = list(sample_dir.glob("*.cif"))
            cif = alt[0] if alt else None
        if cif:
            files["samples"][sample_dir.name] = cif
            # per-sample confidence JSON fallback
            if files["summary"] is None:
                alt_summary = list(sample_dir.glob("*summary_confidences.json"))
                if alt_summary:
                    files["summary"] = alt_summary[0]
            if files["full"] is None:
                alt_full = list(sample_dir.glob("*confidences.json"))
                # exclude the summary file if glob caught both
                alt_full = [f for f in alt_full if "summary" not in f.name]
                if alt_full:
                    files["full"] = alt_full[0]

    return files


def _best_sample_name(ranking_csv: Optional[Path]) -> Optional[str]:
    """Return 'seed-<n>_sample-<m>' for the highest ranking_score row."""
    if ranking_csv is None or not ranking_csv.exists():
        return None
    df = pd.read_csv(ranking_csv)
    if df.empty:
        return None
    best = df.loc[df["ranking_score"].idxmax()]
    return f"seed-{int(best['seed'])}_sample-{int(best['sample'])}"


# ---------------------------------------------------------------------------
# Chain-letter mapping: AF3's *_data.json lists sequences/ligands in the order
# they get assigned chain letters (A, B, C, D, ...) in the output mmCIF. This
# resolves "which JSON confidence index corresponds to which mmCIF chain".
# ---------------------------------------------------------------------------

def get_chain_order(job_dir: Path, job_name: str) -> list[dict]:
    """Return input entities in chain order, e.g.
    [{'letter': 'A', 'kind': 'protein', 'id': 'A'}, {'letter': 'C', 'kind': 'ligand', 'ccd': 'GTP'}, ...]
    """
    data_json = job_dir / f"{job_name}_data.json"
    if not data_json.exists():
        raise FileNotFoundError(f"Could not find {data_json} to resolve chain order")

    with open(data_json) as f:
        spec = json.load(f)

    letters = [chr(ord("A") + i) for i in range(len(spec["sequences"]))]
    entities = []
    for letter, entry in zip(letters, spec["sequences"]):
        if "protein" in entry:
            entities.append({"letter": letter, "kind": "protein", "id": entry["protein"].get("id", letter)})
        elif "ligand" in entry:
            entities.append({
                "letter": letter, "kind": "ligand",
                "ccd": entry["ligand"].get("ccdCodes", ["?"])[0],
            })
        else:
            entities.append({"letter": letter, "kind": "unknown"})
    return entities


# ---------------------------------------------------------------------------
# Structural checks (MDAnalysis) -- equivalents of the ChimeraX commands
# ---------------------------------------------------------------------------

def get_interface_residues(cif_path: Path, chain_a: str, chain_b: str, cutoff: float = 5.0):
    """Equivalent of: select /A & (/B :< cutoff)  -->  returns (resids_a, resids_b)."""
    u = mda.Universe(str(cif_path))
    sel_a = u.select_atoms(f"segid {chain_a} and around {cutoff} segid {chain_b}")
    sel_b = u.select_atoms(f"segid {chain_b} and around {cutoff} segid {chain_a}")
    resids_a = sorted(set(sel_a.resids.tolist()))
    resids_b = sorted(set(sel_b.resids.tolist()))
    return resids_a, resids_b


def get_region_overlap(cif_path: Path, chain_a: str, chain_b: str, ligand_chain: str,
                        interface_cutoff: float = 5.0, ligand_cutoff: float = 10.0):
    """Equivalent of: select /A & (/B :< 5) & (/<ligand> :< 10)
    Residues on chain_a that are BOTH at the chain_a/chain_b interface AND near
    the ligand (e.g. GTP) -- the nucleotide-state-sensitivity check.
    """
    u = mda.Universe(str(cif_path))
    sel = u.select_atoms(
        f"segid {chain_a} and around {interface_cutoff} segid {chain_b} "
        f"and around {ligand_cutoff} segid {ligand_chain}"
    )
    return sorted(set(sel.resids.tolist()))


def compute_buried_area(cif_path: Path, chain_a: str, chain_b: str) -> Optional[float]:
    """Equivalent of ChimeraX `measure buriedarea`. Requires freesasa.
    Writes temp single-chain and complex PDBs and computes:
        buried = (SASA_A + SASA_B - SASA_complex) / 2
    """
    if freesasa is None:
        raise ImportError("freesasa is required for buried-area calc: pip install freesasa --break-system-packages")

    u = mda.Universe(str(cif_path))
    tmpdir = Path(tempfile.mkdtemp(prefix="pocket_verify_"))
    try:
        paths = {}
        for label, sel_str in [
            ("complex", f"segid {chain_a} or segid {chain_b}"),
            ("a", f"segid {chain_a}"),
            ("b", f"segid {chain_b}"),
        ]:
            sel = u.select_atoms(sel_str)
            out = tmpdir / f"{label}.pdb"
            sel.write(str(out))
            paths[label] = out

        areas = {}
        for label, p in paths.items():
            struct = freesasa.Structure(str(p))
            result = freesasa.calc(struct)
            areas[label] = result.totalArea()

        buried = (areas["a"] + areas["b"] - areas["complex"]) / 2.0
        return buried
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Confidence JSON parsing
# ---------------------------------------------------------------------------

def _load_json(path: Optional[Path]) -> dict:
    if path is None or not Path(path).exists():
        return {}
    with open(path) as f:
        return json.load(f)


def get_pairwise_iptm(summary_conf: dict, idx_a: int, idx_b: int) -> Optional[float]:
    """chain_pair_iptm is a square matrix indexed in the same order as the
    input sequences (0-indexed). idx_a/idx_b come from get_chain_order().
    """
    mat = summary_conf.get("chain_pair_iptm")
    if mat is None:
        return None
    try:
        return mat[idx_a][idx_b]
    except (IndexError, TypeError):
        return None


# ---------------------------------------------------------------------------
# Verdict rubric (guide.md step 7, made concrete)
# ---------------------------------------------------------------------------

def call_verdict(pairwise_iptm: Optional[float], buried_area: Optional[float],
                  has_clash: Optional[bool], iptm_std: Optional[float]) -> str:
    if pairwise_iptm is None:
        return "unscored (missing ipTM)"

    if has_clash:
        return "unstable (steric clash flagged by AF3)"

    if buried_area is not None and buried_area < MIN_BURIED_AREA_A2:
        return "weak (interface too small to trust, <400 A^2)"

    if pairwise_iptm >= IPTM_STRONG:
        base = "strong"
    elif pairwise_iptm >= IPTM_WEAK:
        base = "weak"
    else:
        return "unstable (pairwise ipTM < 0.6)"

    if iptm_std is not None and iptm_std > CONSISTENCY_IPTM_STD_MAX:
        return f"{base}, but inconsistent across samples (std={iptm_std:.2f}) -- re-check"

    return base


# ---------------------------------------------------------------------------
# Top-level: score a single job
# ---------------------------------------------------------------------------

def score_job(job_dir: str | Path, chain_a: str, chain_b: str,
              interface_cutoff: float = 5.0, compute_area: bool = True) -> JobScore:
    job_dir = Path(job_dir)
    job_name = job_dir.name

    entities = get_chain_order(job_dir, job_name)
    letter_to_idx = {e["letter"]: i for i, e in enumerate(entities)}
    if chain_a not in letter_to_idx or chain_b not in letter_to_idx:
        raise ValueError(
            f"chain_a/chain_b ({chain_a}, {chain_b}) not found in resolved chain order: {entities}"
        )
    idx_a, idx_b = letter_to_idx[chain_a], letter_to_idx[chain_b]

    files = _find_confidence_files(job_dir, job_name)
    summary_conf = _load_json(files["summary"])
    best_sample_name = _best_sample_name(files["ranking_csv"]) or (
        next(iter(files["samples"]), None)
    )

    result = JobScore(job_name=job_name, chain_a=chain_a, chain_b=chain_b,
                       top_sample_dir=best_sample_name)

    result.global_iptm = summary_conf.get("iptm")
    result.global_ptm = summary_conf.get("ptm")
    result.has_clash = bool(summary_conf.get("has_clash", 0))
    result.fraction_disordered = summary_conf.get("fraction_disordered")
    result.pairwise_iptm = get_pairwise_iptm(summary_conf, idx_a, idx_b)

    # ranking score for the chosen sample
    if files["ranking_csv"] is not None and files["ranking_csv"].exists():
        rdf = pd.read_csv(files["ranking_csv"])
        result.ranked_iptm_values = rdf["ranking_score"].tolist()
        if best_sample_name is not None:
            try:
                seed_str, sample_str = best_sample_name.split("_")
                seed_n, sample_n = int(seed_str.split("-")[1]), int(sample_str.split("-")[1])
                row = rdf[(rdf["seed"] == seed_n) & (rdf["sample"] == sample_n)]
                if not row.empty:
                    result.ranking_score = float(row.iloc[0]["ranking_score"])
            except Exception:
                pass

    iptm_std = float(np.std(result.ranked_iptm_values)) if result.ranked_iptm_values else None

    # structural checks on the best sample's model.cif
    if best_sample_name and best_sample_name in files["samples"]:
        cif_path = files["samples"][best_sample_name]
        try:
            resids_a, resids_b = get_interface_residues(cif_path, chain_a, chain_b, interface_cutoff)
            result.interface_residues_a = resids_a
            result.interface_residues_b = resids_b
            result.n_interface_residues_a = len(resids_a)
            result.n_interface_residues_b = len(resids_b)
        except Exception as e:
            result.notes += f"interface calc failed: {e}; "

        if compute_area:
            try:
                result.buried_area_a2 = compute_buried_area(cif_path, chain_a, chain_b)
            except Exception as e:
                result.notes += f"buried area calc failed: {e}; "
    else:
        result.notes += "no model.cif found for top-ranked sample; "

    result.verdict = call_verdict(result.pairwise_iptm, result.buried_area_a2,
                                   result.has_clash, iptm_std)
    return result


# ---------------------------------------------------------------------------
# Batch scoring across a directory of many job folders (the actual screen)
# ---------------------------------------------------------------------------

def batch_score(base_dir: str | Path, chain_a: str, chain_b: str,
                 job_filter: Optional[list[str]] = None) -> pd.DataFrame:
    """Score every AF3 job subfolder under base_dir. job_filter, if given,
    restricts to job directory names containing any of the given substrings
    (handy for e.g. only scoring FAM129B jobs, or only Tier-1 controls).
    """
    base_dir = Path(base_dir)
    rows = []
    for job_dir in sorted(p for p in base_dir.iterdir() if p.is_dir()):
        if job_filter and not any(f in job_dir.name for f in job_filter):
            continue
        try:
            result = score_job(job_dir, chain_a=chain_a, chain_b=chain_b)
            rows.append(result.as_row())
        except Exception as e:
            rows.append({"job_name": job_dir.name, "verdict": "ERROR", "notes": str(e)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# GTP vs GDP interface comparison -- the nucleotide-state-selectivity check
# ---------------------------------------------------------------------------

def compare_nucleotide_states(gtp_job_dir: str | Path, gdp_job_dir: str | Path,
                               chain_a: str, chain_b: str, interface_cutoff: float = 5.0) -> dict:
    """Compare the RAB-side interface residue set between GTP-bound and
    GDP-bound jobs. A real effector should show interface residues that shrink
    or disappear (especially any switch-region residues) in the GDP state.
    """
    gtp_score = score_job(gtp_job_dir, chain_a=chain_a, chain_b=chain_b, compute_area=False)
    gdp_score = score_job(gdp_job_dir, chain_a=chain_a, chain_b=chain_b, compute_area=False)

    gtp_set = set(gtp_score.interface_residues_a)
    gdp_set = set(gdp_score.interface_residues_a)

    return {
        "gtp_pairwise_iptm": gtp_score.pairwise_iptm,
        "gdp_pairwise_iptm": gdp_score.pairwise_iptm,
        "gtp_only_residues": sorted(gtp_set - gdp_set),
        "gdp_only_residues": sorted(gdp_set - gtp_set),
        "shared_residues": sorted(gtp_set & gdp_set),
        "gtp_interface_size": len(gtp_set),
        "gdp_interface_size": len(gdp_set),
        "interpretation": (
            "Consistent with nucleotide-state-selective binding (interface shrinks in GDP state)"
            if len(gdp_set) < len(gtp_set) and (gtp_score.pairwise_iptm or 0) > (gdp_score.pairwise_iptm or 0)
            else "No clear nucleotide-state selectivity signal -- re-examine before trusting this as a real effector interaction"
        ),
    }


# ---------------------------------------------------------------------------
# CLI entrypoint -- for SLURM job scripts / quick command-line use
# ---------------------------------------------------------------------------

def _cli():
    import argparse

    parser = argparse.ArgumentParser(description="Score AF3 cofold interface(s) for binding-site verification.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_one = sub.add_parser("score", help="Score a single AF3 job directory")
    p_one.add_argument("job_dir")
    p_one.add_argument("--chain-a", required=True)
    p_one.add_argument("--chain-b", required=True)

    p_batch = sub.add_parser("batch", help="Score every job folder under a directory")
    p_batch.add_argument("base_dir")
    p_batch.add_argument("--chain-a", required=True)
    p_batch.add_argument("--chain-b", required=True)
    p_batch.add_argument("--filter", nargs="*", default=None,
                          help="Only include job dirs whose name contains any of these substrings")
    p_batch.add_argument("--out-csv", default="results/scores_summary.csv")

    p_diff = sub.add_parser("compare-nucleotide", help="Compare GTP vs GDP interfaces")
    p_diff.add_argument("gtp_job_dir")
    p_diff.add_argument("gdp_job_dir")
    p_diff.add_argument("--chain-a", required=True)
    p_diff.add_argument("--chain-b", required=True)

    args = parser.parse_args()

    if args.cmd == "score":
        result = score_job(args.job_dir, chain_a=args.chain_a, chain_b=args.chain_b)
        for k, v in result.as_row().items():
            print(f"{k}: {v}")

    elif args.cmd == "batch":
        df = batch_score(args.base_dir, chain_a=args.chain_a, chain_b=args.chain_b, job_filter=args.filter)
        os.makedirs(os.path.dirname(args.out_csv) or ".", exist_ok=True)
        df.to_csv(args.out_csv, index=False)
        print(df[["job_name", "pairwise_iptm", "buried_area_a2", "verdict"]].to_string(index=False))
        print(f"\nWrote full results to {args.out_csv}")

    elif args.cmd == "compare-nucleotide":
        diff = compare_nucleotide_states(args.gtp_job_dir, args.gdp_job_dir,
                                          chain_a=args.chain_a, chain_b=args.chain_b)
        for k, v in diff.items():
            print(f"{k}: {v}")


if __name__ == "__main__":
    _cli()