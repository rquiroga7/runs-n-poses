#!/usr/bin/env python
"""
Generic analysis script for any docking method.

Reads docked outputs from a directory tree
``<output_dir>/<system_id>/<system_id>_<ligand_chain>/`` and produces
``predictions/<method_name>.csv`` + ``posebusters_results/<method_name>.csv``
by running ``ost compare-ligand-structures`` and ``PoseBusters`` for each
ligand pose, mirroring ``04_analyze_vinardo.py``.

The script supports three method families (selected with ``--family``):

* ``vina``      – out.pdbqt, log.txt with VINA RESULT / Affinity lines
* ``vinardo``   – *.pdbt, log.csv with ``Ligand,nConfs,vinardo score``
* ``adgpu``     – output.dlg + output-best.pdbqt, parse first DLG cluster

Usage example:
    python 04_analyze_docking.py \\
        --family vina --method qvina_w \\
        --output-dir /home/rquiroga/Datasets/runs-n-poses-datasets/qvina_w_symmetry \\
        --analysis-dir /home/rquiroga/github/runs-n-poses/examples/analysis/qvina_w \\
        --predictions-csv /home/rquiroga/Datasets/runs-n-poses-datasets/predictions/qvina_w.csv \\
        --posebusters-csv /home/rquiroga/Datasets/runs-n-poses-datasets/posebusters_results/qvina_w.csv \\
        --system-list scripts/systems_for_symmetry_docking.txt
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import inchi

GT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth")
ANNOTATIONS = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/annotations.csv")
OST = "/home/rquiroga/anaconda3/envs/runs_n_poses/bin/ost"
PB_VENV_PYTHON = "/home/rquiroga/Datasets/Posebusters/venv/bin/python"
RDKIT_PYTHON = "/home/rquiroga/anaconda3/envs/RDKIT/bin/python"
SCRIPT_DIR = Path(__file__).parent
FIXER_SCRIPT = SCRIPT_DIR / "fix_pred_copy_gt_valence.py"
COORD_FIXER_SCRIPT = SCRIPT_DIR / "fix_pred_by_coord_match.py"

CANONICAL_PB_COLUMNS = [
    "molecule", "position", "mol_pred_loaded", "mol_true_loaded", "mol_cond_loaded",
    "sanitization", "inchi_convertible", "all_atoms_connected",
    "molecular_formula", "molecular_bonds",
    "double_bond_stereochemistry", "tetrahedral_chirality",
    "bond_lengths", "bond_angles", "internal_steric_clash",
    "aromatic_ring_flatness", "non-aromatic_ring_non-flatness",
    "double_bond_flatness", "internal_energy",
    "protein-ligand_maximum_distance", "minimum_distance_to_protein",
    "minimum_distance_to_organic_cofactors",
    "minimum_distance_to_inorganic_cofactors", "minimum_distance_to_waters",
    "volume_overlap_with_protein", "volume_overlap_with_organic_cofactors",
    "volume_overlap_with_inorganic_cofactors", "volume_overlap_with_waters",
    "system_id", "seed", "sample", "ligand_chain", "method",
]
PB_CHECK_COLUMNS = {
    "sanitization", "inchi_convertible", "all_atoms_connected",
    "bond_lengths", "bond_angles", "internal_steric_clash",
    "aromatic_ring_flatness", "double_bond_flatness",
    "internal_energy", "protein-ligand_maximum_distance",
    "minimum_distance_to_protein",
}


# ---------------------------------------------------------------------------
# CSV utilities
# ---------------------------------------------------------------------------
def _append_unique_csv_rows(csv_path: Path, df: pd.DataFrame) -> int:
    if df.empty:
        return 0

    file_has_header = csv_path.exists() and csv_path.stat().st_size > 0
    if file_has_header:
        with csv_path.open(newline="") as f:
            reader = csv.reader(f)
            try:
                header = next(reader)
            except StopIteration:
                header = list(df.columns)
                existing_rows = set()
                file_has_header = False
            else:
                existing_rows = {tuple(row) for row in reader}
    else:
        header = list(df.columns)
        existing_rows = set()

    df_to_write = df.reindex(columns=header)
    new_rows = []
    for row in df_to_write.itertuples(index=False, name=None):
        row_tuple = tuple("" if pd.isna(value) else str(value) for value in row)
        if row_tuple in existing_rows:
            continue
        existing_rows.add(row_tuple)
        new_rows.append(row_tuple)

    if not new_rows:
        return 0
    with csv_path.open("a", newline="") as f:
        writer = csv.writer(f)
        if not file_has_header:
            writer.writerow(header)
        writer.writerows(new_rows)
    return len(new_rows)


def _coerce_pb_bool(value) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"true", "t", "1", "yes", "y"}:
            return True
        if lowered in {"false", "f", "0", "no", "n", ""}:
            return False
    return bool(value)


# ---------------------------------------------------------------------------
# Pose parsing helpers
# ---------------------------------------------------------------------------
def _convert_pdbt_to_sdf(pdbt_file: str, sdf_file: str) -> bool:
    try:
        r = subprocess.run(
            ["obabel-25-07", pdbt_file, "-O", sdf_file, "-h"],
            capture_output=True, text=True, timeout=30,
        )
        return r.returncode == 0 and os.path.exists(sdf_file) and os.path.getsize(sdf_file) > 100
    except Exception:
        return False


def _parse_vina_score_from_log(log_file: Path) -> float:
    """Best-effort score extraction from a Vina-family log.txt."""
    if not log_file.exists():
        return None
    try:
        for ln in reversed(log_file.read_text().strip().split("\n")[-20:]):
            if "VINA RESULT" in ln or "Affinity" in ln:
                parts = ln.replace(",", " ").split()
                for p in parts:
                    try:
                        return float(p)
                    except Exception:
                        continue
    except Exception:
        pass
    return None


def _parse_vinardo_score_from_log(log_file: Path) -> float:
    if not log_file.exists():
        return None
    try:
        with open(log_file) as f:
            lines = f.read().strip().split("\n")
        if len(lines) < 2:
            return None
        parts = lines[1].split(",")
        if len(parts) >= 3:
            return float(parts[2])
    except Exception:
        pass
    return None


def _parse_adgpu_score(dlg_file: Path) -> float:
    """Parse the best energy from an AutoDock DLG file."""
    if not dlg_file.exists():
        return None
    try:
        with open(dlg_file) as f:
            text = f.read()
        # Energy line format: "    1 |        -7.123 |    ...". We take the first one.
        for ln in text.splitlines():
            if "lowest" in ln.lower() and "rms" not in ln.lower():
                continue
            if "|" in ln and "kcal/mol" in ln:
                parts = [p for p in ln.split("|") if p.strip()]
                if parts:
                    try:
                        return float(parts[0].strip())
                    except Exception:
                        continue
    except Exception:
        pass
    return None


# ---------------------------------------------------------------------------
# PoseBusters + RDKit fixer helpers (reused from 04_analyze_vinardo.py)
# ---------------------------------------------------------------------------
def _run_rdkit_fixer(python_exe: str, script: Path, pred_path: str, gt_path: str) -> str | None:
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
    tmp.close()
    try:
        cmd = [python_exe, str(script), pred_path, gt_path, tmp.name]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if r.returncode == 0 and os.path.exists(tmp.name) and os.path.getsize(tmp.name) > 0:
            return tmp.name
    except Exception:
        pass
    if os.path.exists(tmp.name):
        os.remove(tmp.name)
    return None


def _is_inchi_convertible_sdf(sdf_path: str) -> bool:
    mol = Chem.MolFromMolFile(sdf_path, sanitize=False, removeHs=False)
    if mol is None:
        return False
    try:
        Chem.SanitizeMol(mol)
        inchi.MolToInchiKey(mol)
        return True
    except Exception:
        return False


def build_posebusters_pred_input(docked_pose: str, docked_sdf: str, gt_sdf: str) -> tuple[str, list[str]]:
    cleanup: list[str] = []

    repaired = _run_rdkit_fixer(RDKIT_PYTHON, FIXER_SCRIPT, docked_sdf, gt_sdf)
    if repaired:
        cleanup.append(repaired)
        if _is_inchi_convertible_sdf(repaired):
            return repaired, cleanup

    pdbt_sdf = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
    pdbt_sdf.close()
    if _convert_pdbt_to_sdf(docked_pose, pdbt_sdf.name):
        cleanup.append(pdbt_sdf.name)
        repaired = _run_rdkit_fixer(RDKIT_PYTHON, COORD_FIXER_SCRIPT, pdbt_sdf.name, gt_sdf)
        if repaired:
            cleanup.append(repaired)
            if _is_inchi_convertible_sdf(repaired):
                return repaired, cleanup
            fallback = _run_rdkit_fixer(RDKIT_PYTHON, FIXER_SCRIPT, repaired, gt_sdf)
            if fallback:
                cleanup.append(fallback)
                return fallback, cleanup
    elif os.path.exists(pdbt_sdf.name):
        os.remove(pdbt_sdf.name)

    repaired = _run_rdkit_fixer(RDKIT_PYTHON, COORD_FIXER_SCRIPT, docked_pose, gt_sdf)
    if repaired:
        cleanup.append(repaired)
        if _is_inchi_convertible_sdf(repaired):
            return repaired, cleanup
        fallback = _run_rdkit_fixer(RDKIT_PYTHON, FIXER_SCRIPT, repaired, gt_sdf)
        if fallback:
            cleanup.append(fallback)
            return fallback, cleanup

    return docked_sdf, cleanup


def run_ost_comparison(receptor_cif: str, docked_sdf: str, gt_sdf: str, output_json: str) -> bool:
    cmd = [
        OST, "compare-ligand-structures",
        "-m", receptor_cif, "-ml", docked_sdf,
        "-r", receptor_cif, "-rl", gt_sdf,
        "-o", output_json, "--lddt-pli", "--rmsd", "--lddt-pli-amc",
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        return r.returncode == 0 and os.path.exists(output_json)
    except Exception:
        return False


def extract_metrics_from_ost(json_file: str, system_id: str, ligand_chain: str,
                              annotations_row: dict, ranking_score, method: str) -> dict | None:
    try:
        with open(json_file) as f:
            result = json.load(f)
    except Exception:
        return None
    if result.get("status") != "SUCCESS":
        return None

    metrics = {
        "target": system_id,
        "method": method,
        "seed": 1,
        "sample": 1,
        "ranking_score": ranking_score,
        "ligand_instance_chain": ligand_chain,
        "ligand_is_proper": annotations_row.get("ligand_is_proper", True),
        "model_ligand_ccd_code": annotations_row.get("ligand_ccd_code", ""),
        "model_ligand_smiles": annotations_row.get("ligand_smiles", ""),
        "ligand_ccd_code": annotations_row.get("ligand_ccd_code", ""),
    }
    if "lddt_pli" in result and "assigned_scores" in result["lddt_pli"]:
        for item in result["lddt_pli"]["assigned_scores"]:
            metrics["lddt_pli"] = item.get("score")
            break
    if "rmsd" in result and "assigned_scores" in result["rmsd"]:
        for item in result["rmsd"]["assigned_scores"]:
            metrics["rmsd"] = item.get("score")
            metrics["lddt_lp"] = item.get("lddt_lp")
            metrics["bb_rmsd"] = item.get("bb_rmsd")
            break
    for col in [
        "prot_lig_chain_iptm_average_lddt_pli", "prot_lig_chain_iptm_min_lddt_pli",
        "prot_lig_chain_iptm_max_lddt_pli", "lig_prot_chain_iptm_average_lddt_pli",
        "lig_prot_chain_iptm_min_lddt_pli", "lig_prot_chain_iptm_max_lddt_pli",
        "prot_lig_chain_iptm_average_rmsd", "prot_lig_chain_iptm_min_rmsd",
        "prot_lig_chain_iptm_max_rmsd", "lig_prot_chain_iptm_average_rmsd",
        "lig_prot_chain_iptm_min_rmsd", "lig_prot_chain_iptm_max_rmsd",
    ]:
        metrics[col] = None
    for col in [
        "pred_pocket_tp", "pred_pocket_fp", "pred_pocket_fn",
        "pred_pocket_precision", "pred_pocket_recall", "pred_pocket_f1",
    ]:
        metrics[col] = None
    return metrics


def run_posebusters(model_sdf: str, gt_sdf: str, receptor_cif: str) -> dict:
    try:
        receptor_pdb = receptor_cif.rsplit(".", 1)[0] + ".pdb"
        subprocess.run(
            ["obabel-25-07", receptor_cif, "-O", receptor_pdb],
            capture_output=True, text=True, timeout=30,
        )
        cmd = [PB_VENV_PYTHON, "-c", f"""
import json
from posebusters import PoseBusters
bust = PoseBusters(config="redock")
results = bust.bust(mol_pred="{model_sdf}", mol_true="{gt_sdf}", mol_cond="{receptor_pdb}")
checks = [
    'sanitization', 'inchi_convertible', 'all_atoms_connected',
    'bond_lengths', 'bond_angles', 'internal_steric_clash',
    'aromatic_ring_flatness', 'double_bond_flatness',
    'internal_energy', 'protein-ligand_maximum_distance',
    'minimum_distance_to_protein',
]
available = [c for c in checks if c in results.columns]
results['pb_success'] = results[available].all(axis=1)
print(results.to_json(orient="records"))
"""]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if r.returncode == 0 and r.stdout.strip():
            results_list = json.loads(r.stdout.strip())
            if results_list:
                return results_list[0]
    except Exception:
        pass
    return {"pb_success": False}


# ---------------------------------------------------------------------------
# Family-specific per-ligand analysis
# ---------------------------------------------------------------------------
def _find_docked_file(lig_dir: Path, family: str) -> Path | None:
    if family == "vina":
        for name in ("out.pdbqt", "out.pdb"):
            p = lig_dir / name
            if p.exists() and p.stat().st_size > 0:
                return p
        return None
    if family == "vinardo":
        candidates = sorted(lig_dir.glob("*.pdbt"))
        for c in candidates:
            if c.stat().st_size > 100:
                return c
        return None
    if family == "adgpu":
        # Prefer the best pose PDBQT; fall back to the DLG.
        for name in ("output-best.pdbqt", "best.pdbqt"):
            p = lig_dir / name
            if p.exists() and p.stat().st_size > 0:
                return p
        p = lig_dir / "output.dlg"
        if p.exists() and p.stat().st_size > 0:
            return p
        return None
    return None


def _find_dlg_for_score(lig_dir: Path) -> Path | None:
    p = lig_dir / "output.dlg"
    return p if p.exists() and p.stat().st_size > 0 else None


def analyze_system(system_id: str, annotations: pd.DataFrame, args) -> tuple[list, list]:
    family = args.family
    method = args.method
    output_dir = Path(args.output_dir)
    analysis_dir = Path(args.analysis_dir)

    system_annot = annotations[annotations["system_id"] == system_id]
    if system_annot.empty:
        return [], []

    results: list = []
    pb_rows: list = []
    system_out = output_dir / system_id
    if not system_out.exists():
        return [], []

    receptor_cif = GT_DIR / system_id / "receptor.cif"
    if not receptor_cif.exists():
        return [], []

    # Build list of (chain, dir) entries — handle both nested and flat output
    subdirs = [d for d in sorted(system_out.iterdir()) if d.is_dir()]
    if subdirs:
        # Nested: <sys_id>/<sys_id>_<chain>/
        lig_entries = []
        for d in subdirs:
            parts = d.name.split("_")
            chain = parts[-1] if len(parts) >= 2 else d.name
            lig_entries.append((chain, d))
    else:
        # Flat: <sys_id>/ files (adgpu style) — get chain from system_annot
        chain = system_annot.iloc[0]["ligand_instance_chain"]
        lig_entries = [(chain, system_out)]

    def _analyze_one(lig_dir: Path, chain: str, method_override: str | None = None) -> tuple[dict | None, dict | None]:
        nonlocal system_id, system_annot, receptor_cif, analysis_dir
        cur_method = method_override or method

        score = None
        if family == "vina":
            score = _parse_vina_score_from_log(lig_dir / "log.txt")
        elif family == "vinardo":
            score = _parse_vinardo_score_from_log(lig_dir / "log.csv")
        elif family == "adgpu":
            dlg = _find_dlg_for_score(lig_dir)
            score = _parse_adgpu_score(dlg) if dlg else None

        docked_file = _find_docked_file(lig_dir, family)
        if docked_file is None:
            return None, None

        sdf_suffix = f"_{cur_method}" if cur_method != method else ""
        docked_sdf = lig_dir / f"{chain}_dock{sdf_suffix}.sdf"
        if not docked_sdf.exists() or docked_sdf.stat().st_size < 100:
            if not _convert_pdbt_to_sdf(str(docked_file), str(docked_sdf)):
                return None, None

        gt_sdf = GT_DIR / system_id / "ligand_files" / f"{chain}.sdf"
        if not gt_sdf.exists():
            return None, None

        analysis_dir.mkdir(parents=True, exist_ok=True)
        json_suffix = f"_{cur_method}" if cur_method != method else ""
        json_file = analysis_dir / f"{system_id}_{chain}{json_suffix}.json"
        metrics = None
        if not json_file.exists():
            if not run_ost_comparison(str(receptor_cif), str(docked_sdf), str(gt_sdf), str(json_file)):
                if args.skip_posebusters:
                    return None, None
        if json_file.exists():
            annot_row = system_annot[system_annot["ligand_instance_chain"] == chain]
            if not annot_row.empty:
                metrics = extract_metrics_from_ost(
                    str(json_file), system_id, chain,
                    annot_row.iloc[0].to_dict(), score, cur_method,
                )

        pb_row = None
        if not args.skip_posebusters:
            pb_pred, tmp_files = build_posebusters_pred_input(
                str(docked_file), str(docked_sdf), str(gt_sdf),
            )
            pb_results = run_posebusters(pb_pred, str(gt_sdf), str(receptor_cif))
            try:
                for tmp_path in tmp_files:
                    if tmp_path and os.path.exists(tmp_path):
                        os.remove(tmp_path)
            except Exception:
                pass

            if metrics is not None:
                metrics["pb_success"] = 1.0 if pb_results.get("pb_success", False) else 0.0
                rt = None
                rt_file = lig_dir / "runtime.json"
                if rt_file.exists():
                    try:
                        with open(rt_file) as f:
                            rt = json.load(f).get("runtime_seconds")
                    except Exception:
                        pass
                if rt is not None:
                    metrics["runtime_seconds"] = rt

            pb_row = dict(pb_results) if isinstance(pb_results, dict) else {}
            pb_row["pb_success"] = 1.0 if pb_results.get("pb_success", False) else 0.0
            pb_row["system_id"] = system_id
            pb_row.setdefault("seed", 1)
            pb_row.setdefault("sample", 1)
            pb_row.setdefault("ligand_chain", chain)
            pb_row.setdefault("method", cur_method)
            rt = None
            rt_file = lig_dir / "runtime.json"
            if rt_file.exists():
                try:
                    with open(rt_file) as f:
                        rt = json.load(f).get("runtime_seconds")
                except Exception:
                    pass
            if rt is not None:
                pb_row["runtime_seconds"] = rt
        else:
            if metrics is not None:
                metrics["pb_success"] = -1.0

        return metrics, pb_row

    for chain, lig_dir in lig_entries:
        metrics, pb_row = _analyze_one(lig_dir, chain)
        if metrics is not None:
            results.append(metrics)
            if pb_row is not None:
                pb_rows.append(pb_row)

        # For adgpu, also analyze the most-populated cluster pose (output-mostpop.pdbqt)
        if family == "adgpu":
            mostpop_path = lig_dir / "output-mostpop.pdbqt"
            if mostpop_path.exists() and mostpop_path.stat().st_size > 0:
                # Temporarily make the analyze function find output-mostpop.pdbqt
                orig_find = _find_docked_file
                def _find_mostpop(l, f):
                    p = l / "output-mostpop.pdbqt"
                    if p.exists() and p.stat().st_size > 0:
                        return p
                    return orig_find(l, f)
                import types
                saved_find = _find_docked_file
                try:
                    # Monkey-patch for mostpop
                    import sys
                    mod = sys.modules[__name__]
                    mod._find_docked_file = _find_mostpop
                    metrics_m, pb_row_m = _analyze_one(lig_dir, chain, method_override=f"{method}_mostpop")
                    if metrics_m is not None:
                        results.append(metrics_m)
                        if pb_row_m is not None:
                            pb_rows.append(pb_row_m)
                finally:
                    mod._find_docked_file = saved_find

    return results, pb_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family", choices=["vina", "vinardo", "adgpu"], required=True)
    parser.add_argument("--method", required=True, help="method name written to predictions CSV")
    parser.add_argument("--output-dir", required=True, help="per-system docking output directory")
    parser.add_argument("--analysis-dir", required=True, help="where to cache OST JSON files")
    parser.add_argument("--predictions-csv", required=True)
    parser.add_argument("--posebusters-csv", required=True)
    parser.add_argument("--system-list", required=True)
    parser.add_argument("--start-index", type=int, default=0)
    parser.add_argument("--max-systems", type=int, default=None)
    parser.add_argument("--skip-posebusters", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    with open(args.system_list) as f:
        system_ids = [l.strip() for l in f if l.strip()]
    system_ids = system_ids[args.start_index:]
    if args.max_systems:
        system_ids = system_ids[:args.max_systems]

    annotations = pd.read_csv(ANNOTATIONS)
    run_pb = not args.skip_posebusters

    pred_csv = Path(args.predictions_csv)
    pb_csv = Path(args.posebusters_csv)
    pred_csv.parent.mkdir(parents=True, exist_ok=True)
    pb_csv.parent.mkdir(parents=True, exist_ok=True)

    processed_pred = set()
    processed_pb = set()
    if args.resume:
        if pred_csv.exists():
            try:
                df = pd.read_csv(pred_csv, low_memory=False)
                if "target" in df.columns:
                    processed_pred = set(df["target"].astype(str).unique())
            except Exception:
                processed_pred = set()
        if pb_csv.exists():
            try:
                df = pd.read_csv(pb_csv, low_memory=False)
                if "system_id" in df.columns:
                    processed_pb = set(df["system_id"].astype(str).unique())
            except Exception:
                processed_pb = set()

    total_metrics = 0
    n_workers = max(1, (os.cpu_count() or 4))
    from concurrent.futures import ThreadPoolExecutor, as_completed
    pool = ThreadPoolExecutor(max_workers=n_workers)

    def _worker(sid):
        if args.resume and sid in processed_pred and sid in processed_pb:
            return sid, [], []
        return sid, *analyze_system(sid, annotations, args)

    futures = {pool.submit(_worker, sid): sid for sid in system_ids}
    for fut in tqdm(as_completed(futures), total=len(system_ids), desc=f"Analyzing {args.method}"):
        try:
            sid, metrics, pb_rows = fut.result()
        except Exception:
            continue
        if not metrics and not pb_rows:
            continue
        if metrics:
            df = pd.DataFrame(metrics)
            appended = _append_unique_csv_rows(pred_csv, df)
            total_metrics += appended
        if pb_rows and run_pb:
            df_pb = pd.DataFrame(pb_rows)
            for c, default in [("system_id", sid), ("seed", 1), ("sample", 1),
                               ("ligand_chain", None), ("method", args.method)]:
                if c not in df_pb.columns:
                    df_pb[c] = default
            for col in CANONICAL_PB_COLUMNS:
                if col not in df_pb.columns:
                    df_pb[col] = pd.NA
            df_pb = df_pb[CANONICAL_PB_COLUMNS]
            for chk in PB_CHECK_COLUMNS:
                if chk in df_pb.columns:
                    df_pb[chk] = df_pb[chk].map(_coerce_pb_bool)
            _append_unique_csv_rows(pb_csv, df_pb)
    pool.shutdown()

    if total_metrics:
        print(f"\nResults saved to: {pred_csv}")
        print(f"  Total predictions: {total_metrics}")


if __name__ == "__main__":
    main()

