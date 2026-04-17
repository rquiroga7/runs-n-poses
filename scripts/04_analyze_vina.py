#!/usr/bin/env python
"""
Step 4 (VINA): Analyze Vina docking results using ost (same as original pipeline).

This is a copy of `04_analyze_vinardo.py` adapted to read Vina outputs and write
predictions/posebusters outputs for the `vina` method.
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm
import tempfile

GT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth")

OUTPUT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/vina_outputs")
ANALYSIS_DIR = Path("/home/rquiroga/github/runs-n-poses/examples/analysis/vina")
ANNOTATIONS = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/annotations.csv")
OUTPUT_CSV = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/predictions/vina.csv")

POSEBUSTERS_CSV = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/posebusters_results/vina.csv")

# Keep other constants same as vinardo analysis
RDKIT_PYTHON = "/home/rquiroga/anaconda3/envs/runs_n_poses/bin/python"
FIXER_SCRIPT = Path(__file__).parent / "fix_pred_copy_gt_valence.py"
CANONICAL_PB_COLUMNS = [
    "molecule",
    "position",
    "mol_pred_loaded",
    "mol_true_loaded",
    "mol_cond_loaded",
    "sanitization",
    "inchi_convertible",
    "all_atoms_connected",
    "molecular_formula",
    "molecular_bonds",
    "double_bond_stereochemistry",
    "tetrahedral_chirality",
    "bond_lengths",
    "bond_angles",
    "internal_steric_clash",
    "aromatic_ring_flatness",
    "non-aromatic_ring_non-flatness",
    "double_bond_flatness",
    "internal_energy",
    "protein-ligand_maximum_distance",
    "minimum_distance_to_protein",
    "minimum_distance_to_organic_cofactors",
    "minimum_distance_to_inorganic_cofactors",
    "minimum_distance_to_waters",
    "volume_overlap_with_protein",
    "volume_overlap_with_organic_cofactors",
    "volume_overlap_with_inorganic_cofactors",
    "volume_overlap_with_waters",
    "system_id",
    "seed",
    "sample",
    "ligand_chain",
    "method",
]

PB_CHECK_COLUMNS = {
    "sanitization", "inchi_convertible", "all_atoms_connected",
    "bond_lengths", "bond_angles", "internal_steric_clash",
    "aromatic_ring_flatness", "double_bond_flatness",
    "internal_energy", "protein-ligand_maximum_distance",
    "minimum_distance_to_protein"
}

OST = "/home/rquiroga/anaconda3/envs/runs_n_poses/bin/ost"
PB_VENV_PYTHON = "/home/rquiroga/Datasets/Posebusters/venv/bin/python"


def convert_pdbt_to_sdf(pdbt_file: str, sdf_file: str) -> bool:
    """Convert PDBT/PDBQT/PDB to SDF using obabel."""
    try:
        r = subprocess.run(
            ["obabel-25-07", pdbt_file, "-O", sdf_file, "-h"],
            capture_output=True, text=True, timeout=30
        )
        return r.returncode == 0 and os.path.exists(sdf_file) and os.path.getsize(sdf_file) > 100
    except Exception:
        return False


def run_ost_comparison(receptor_cif: str, docked_sdf: str, gt_sdf: str, output_json: str) -> bool:
    cmd = [
        OST, "compare-ligand-structures",
        "-m", receptor_cif,
        "-ml", docked_sdf,
        "-r", receptor_cif,
        "-rl", gt_sdf,
        "-o", output_json,
        "--lddt-pli", "--rmsd", "--lddt-pli-amc"
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        return r.returncode == 0 and os.path.exists(output_json)
    except Exception:
        return False


def extract_metrics_from_ost(json_file: str, system_id: str, ligand_chain: str,
                              annotations_row: dict, vina_score: float) -> dict:
    try:
        with open(json_file) as f:
            result = json.load(f)
    except Exception:
        return None

    if result.get("status") != "SUCCESS":
        return None

    metrics = {
        "target": system_id,
        "method": "vina",
        "seed": 1,
        "sample": 1,
        "ranking_score": vina_score,
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

    metrics["pred_pocket_tp"] = None
    metrics["pred_pocket_fp"] = None
    metrics["pred_pocket_fn"] = None
    metrics["pred_pocket_precision"] = None
    metrics["pred_pocket_recall"] = None
    metrics["pred_pocket_f1"] = None

    return metrics


def run_posebusters(model_sdf: str, gt_sdf: str, receptor_cif: str) -> dict:
    try:
        receptor_pdb = receptor_cif.rsplit(".", 1)[0] + ".pdb"
        subprocess.run([
            "obabel-25-07", receptor_cif, "-O", receptor_pdb
        ], capture_output=True, text=True, timeout=30)

        cmd = [
            PB_VENV_PYTHON, "-c",
            f"""
import json
from posebusters import PoseBusters
bust = PoseBusters(config=\"redock\")
results = bust.bust(
    mol_pred=\"{model_sdf}\", 
    mol_true=\"{gt_sdf}\",
    mol_cond=\"{receptor_pdb}\",
)
checks = [
    'sanitization', 'inchi_convertible', 'all_atoms_connected',
    'bond_lengths', 'bond_angles', 'internal_steric_clash',
    'aromatic_ring_flatness', 'double_bond_flatness',
    'internal_energy', 'protein-ligand_maximum_distance',
    'minimum_distance_to_protein'
]
available = [c for c in checks if c in results.columns]
results['pb_success'] = results[available].all(axis=1)
print(results.to_json(orient=\"records\"))
"""
        ]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if r.returncode == 0 and r.stdout.strip():
            results_list = json.loads(r.stdout.strip())
            if results_list:
                return results_list[0]
    except Exception:
        pass
    return {"pb_success": False}


# Reuse RDKit sanitizer helpers and analyze_system logic from vinardo script
# For brevity we'll import the sanitizer helper if available by relative import
try:
    from .fix_pred_copy_gt_valence import fix_pred_copy_gt_valence
except Exception:
    pass


# For make_rdkit_sanitized_sdf we will reuse the same implementation pattern
try:
    from rdkit import Chem
except Exception:
    Chem = None


def make_rdkit_sanitized_sdf(src_path: str) -> str:
    if Chem is None:
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        tmp.close()
        try:
            r = subprocess.run([
                "obabel-25-07", src_path, "-O", tmp.name, "-h"
            ], capture_output=True, text=True, timeout=30)
            if r.returncode == 0 and os.path.exists(tmp.name) and os.path.getsize(tmp.name) > 0:
                return tmp.name
        except Exception:
            pass
        if os.path.exists(tmp.name):
            os.remove(tmp.name)
        raise RuntimeError("Neither RDKit available nor OpenBabel round-trip succeeded")

    suppl = Chem.SDMolSupplier(src_path, removeHs=False)
    mols = [m for m in suppl if m is not None]
    if not mols:
        suppl = Chem.SDMolSupplier(src_path, removeHs=True)
        mols = [m for m in suppl if m is not None]
    if not mols:
        raise RuntimeError(f"RDKit could not read molecules from: {src_path}")
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
    tmp.close()
    writer = Chem.SDWriter(tmp.name)
    for m in mols:
        if m is None:
            continue
        try:
            Chem.SanitizeMol(m)
        except Exception:
            try:
                Chem.SanitizeMol(m, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
            except Exception:
                pass
        try:
            m_h = Chem.AddHs(m)
            writer.write(m_h)
        except Exception:
            try:
                writer.write(m)
            except Exception:
                pass
    writer.close()
    if not os.path.exists(tmp.name) or os.path.getsize(tmp.name) == 0:
        if os.path.exists(tmp.name):
            os.remove(tmp.name)
        raise RuntimeError("Failed to write RDKit-sanitized SDF")
    return tmp.name


# The analyze_system implementation is intentionally similar to vinardo's but
# adjusted to handle Vina output naming (out.pdbqt / log.txt) and method name 'vina'.

def analyze_system(system_id: str, annotations: pd.DataFrame, run_posebusters_flag: bool = True,
                   run_ost_flag: bool = True):
    system_annot = annotations[annotations["system_id"] == system_id]
    if system_annot.empty:
        return [], []

    results = []
    pb_rows = []
    system_out = OUTPUT_DIR / system_id
    if not system_out.exists():
        return [], []

    receptor_cif = GT_DIR / system_id / "receptor.cif"
    if not receptor_cif.exists():
        return [], []

    for lig_dir in sorted(system_out.iterdir()):
        if not lig_dir.is_dir():
            continue

        parts = lig_dir.name.split("_")
        if len(parts) < 2:
            continue
        chain = parts[-1]

        # Check if docking succeeded: look for log.txt or out.pdbqt
        log_file = lig_dir / "log.txt"
        out_pdbqt = lig_dir / "out.pdbqt"
        if not (log_file.exists() and log_file.stat().st_size > 0) and not out_pdbqt.exists():
            continue

        # Extract vina score (attempt from log)
        vina_score = None
        try:
            if log_file.exists():
                log_lines = log_file.read_text().strip().split("\n")
                for ln in reversed(log_lines[-20:]):
                    if "Affinity" in ln or "REMARK" in ln:
                        # best-effort parse; not guaranteed
                        parts = ln.replace(',', ' ').split()
                        for p in parts:
                            try:
                                v = float(p)
                                vina_score = v
                                break
                            except Exception:
                                continue
                        if vina_score is not None:
                            break
        except Exception:
            pass

        # Find docked pdbqt/pdb (prefer out.pdbqt)
        docked_file = None
        candidates = list(lig_dir.glob("*.pdbqt")) + list(lig_dir.glob("*.pdb"))
        if candidates:
            docked_file = candidates[0]
        else:
            # nothing to convert
            continue

        # Convert to SDF
        docked_sdf = lig_dir / f"{chain}_dock.sdf"
        if not docked_sdf.exists() or docked_sdf.stat().st_size < 100:
            if not convert_pdbt_to_sdf(str(docked_file), str(docked_sdf)):
                continue

        gt_sdf = GT_DIR / system_id / "ligand_files" / f"{chain}.sdf"
        if not gt_sdf.exists():
            continue

        ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
        json_file = ANALYSIS_DIR / f"{system_id}_{chain}.json"

        metrics = None
        if run_ost_flag:
            if not json_file.exists():
                if not run_ost_comparison(str(receptor_cif), str(docked_sdf), str(gt_sdf), str(json_file)):
                    if not run_posebusters_flag:
                        continue
            if json_file.exists():
                annot_row = system_annot[system_annot["ligand_instance_chain"] == chain]
                if not annot_row.empty:
                    metrics = extract_metrics_from_ost(str(json_file), system_id, chain,
                                                       annot_row.iloc[0].to_dict(), vina_score)
                else:
                    metrics = None

        if run_posebusters_flag:
            tmp_pred = None
            tmp_gt = None
            tmp_fixed_pred = None
            try:
                tmp_pred = make_rdkit_sanitized_sdf(str(docked_sdf))
            except Exception:
                tmp_pred = None
            try:
                tmp_gt = make_rdkit_sanitized_sdf(str(gt_sdf))
            except Exception:
                tmp_gt = None

            try:
                tf = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
                tf.close()
                fixer_cmd = [RDKIT_PYTHON, str(FIXER_SCRIPT), str(docked_sdf), str(gt_sdf), tf.name]
                r = subprocess.run(fixer_cmd, capture_output=True, text=True, timeout=30)
                if r.returncode == 0 and os.path.exists(tf.name) and os.path.getsize(tf.name) > 0:
                    tmp_fixed_pred = tf.name
                else:
                    if os.path.exists(tf.name):
                        os.remove(tf.name)
                    tmp_fixed_pred = None
            except Exception:
                tmp_fixed_pred = None

            use_pred = tmp_fixed_pred if tmp_fixed_pred else (tmp_pred if tmp_pred else str(docked_sdf))
            use_gt = tmp_gt if tmp_gt else str(gt_sdf)
            pb_results = run_posebusters(use_pred, use_gt, str(receptor_cif))
            try:
                if tmp_pred and os.path.exists(tmp_pred):
                    os.remove(tmp_pred)
                if tmp_gt and os.path.exists(tmp_gt):
                    os.remove(tmp_gt)
                if tmp_fixed_pred and os.path.exists(tmp_fixed_pred):
                    os.remove(tmp_fixed_pred)
            except Exception:
                pass

            if metrics is not None:
                metrics["pb_success"] = 1.0 if pb_results.get("pb_success", False) else 0.0

            pb_row = dict(pb_results) if isinstance(pb_results, dict) else {}
            try:
                pb_row["pb_success"] = 1.0 if pb_results.get("pb_success", False) else 0.0
            except Exception:
                pb_row["pb_success"] = 0.0
            pb_row["system_id"] = system_id
            pb_row.setdefault("seed", (metrics.get("seed", 1) if metrics is not None else 1))
            pb_row.setdefault("sample", (metrics.get("sample", 1) if metrics is not None else 1))
            pb_row.setdefault("ligand_chain", chain)
            pb_row.setdefault("method", "vina")
            pb_rows.append(pb_row)
        else:
            if metrics is not None:
                metrics["pb_success"] = -1.0

        if metrics is not None:
            results.append(metrics)

    return results, pb_rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--system-list", type=str, default="scripts/single_ligand_systems.txt")
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
    resume = args.resume

    processed_pred_targets = set()
    processed_pb_targets = set()
    if resume:
        if OUTPUT_CSV.exists():
            try:
                existing_preds = pd.read_csv(OUTPUT_CSV, low_memory=False)
                if 'target' in existing_preds.columns:
                    processed_pred_targets = set(existing_preds['target'].astype(str).unique())
            except Exception:
                processed_pred_targets = set()
        if POSEBUSTERS_CSV.exists():
            try:
                existing_pb = pd.read_csv(POSEBUSTERS_CSV, low_memory=False)
                if 'system_id' in existing_pb.columns:
                    processed_pb_targets = set(existing_pb['system_id'].astype(str).unique())
            except Exception:
                processed_pb_targets = set()

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    POSEBUSTERS_CSV.parent.mkdir(parents=True, exist_ok=True)

    total_metrics = 0
    for sid in tqdm(system_ids, desc="Analyzing (Vina)"):
        pred_done = resume and (sid in processed_pred_targets)
        pb_done = resume and (sid in processed_pb_targets)
        if pred_done and pb_done:
            continue
        run_ost_flag = not pred_done
        run_pb_for_system = run_pb and (not pb_done)

        metrics, pb_rows = analyze_system(sid, annotations,
                                         run_posebusters_flag=run_pb_for_system,
                                         run_ost_flag=run_ost_flag)

        if (not metrics) and (not pb_rows):
            continue

        if metrics and run_ost_flag:
            df_sys = pd.DataFrame(metrics)
            write_header = not OUTPUT_CSV.exists()
            try:
                df_sys.to_csv(OUTPUT_CSV, index=False, header=write_header, mode='a')
            except Exception:
                df_sys.to_csv(OUTPUT_CSV, index=False)
            processed_pred_targets.add(sid)

        if pb_rows and run_pb_for_system:
            df_pb = pd.DataFrame(pb_rows)
            for c, default in [("system_id", sid), ("seed", 1), ("sample", 1), ("ligand_chain", None), ("method", "vina")]:
                if c not in df_pb.columns:
                    df_pb[c] = default
            for col in CANONICAL_PB_COLUMNS:
                if col not in df_pb.columns:
                    df_pb[col] = pd.NA
            df_pb = df_pb[CANONICAL_PB_COLUMNS]
            for chk in PB_CHECK_COLUMNS:
                if chk in df_pb.columns:
                    df_pb[chk] = df_pb[chk].fillna(False)
            write_header_pb = not POSEBUSTERS_CSV.exists()
            try:
                df_pb.to_csv(POSEBUSTERS_CSV, index=False, header=write_header_pb, mode='a')
            except Exception:
                df_pb.to_csv(POSEBUSTERS_CSV, index=False)
            processed_pb_targets.add(sid)

        if metrics and run_ost_flag:
            total_metrics += len(df_sys)

    if total_metrics > 0:
        df_final = pd.read_csv(OUTPUT_CSV)
        print(f"\nResults saved to: {OUTPUT_CSV}")
        print(f"  Total predictions: {len(df_final)}")
        print(f"  Unique systems: {df_final['target'].nunique()}")
    else:
        print("No new metrics extracted.")


if __name__ == "__main__":
    main()
