#!/usr/bin/env python
"""
Step 4: Analyze vinardo docking results using ost (same as original pipeline).

For each successfully docked system:
1. Convert PDBT→SDF for the docked ligand
2. Run ost compare-ligand-structures with -ml flag:
   - Model: ground truth receptor + docked ligand SDF (-ml)
   - Reference: ground truth receptor + ground truth ligand SDF (-rl)
3. Run PoseBusters for validation checks
4. Extract lddt_pli, rmsd, lddt_lp, bb_rmsd metrics
5. Output CSV in same format as other prediction methods

Usage:
    python 04_analyze_vinardo.py [--start-index N] [--max-systems N] [--skip-posebusters]
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
import atexit


OUTPUT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_outputs")
ANALYSIS_DIR = Path("/home/rquiroga/github/runs-n-poses/examples/analysis/vinardock_2vinardo")
ANNOTATIONS = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/annotations.csv")
OUTPUT_CSV = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/predictions/vinardock_2vinardo.csv")


POSEBUSTERS_CSV = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/posebusters_results/vinardock_2vinardo.csv")

# Python executable for RDKit environment (used by fixer script)
RDKIT_PYTHON = "/home/rquiroga/anaconda3/envs/runs_n_poses/bin/python"
# Path to fixer script that copies GT valence/H onto predicted ligands
FIXER_SCRIPT = Path(__file__).parent / "fix_pred_copy_gt_valence.py"
# Canonical PoseBusters (AF3) column order - used to normalize output CSVs
# Exact header copied from runs-n-poses-datasets/posebusters_results/af3.csv
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

# Subset of check columns that are boolean pass/fail in PoseBusters
PB_CHECK_COLUMNS = {
    "sanitization", "inchi_convertible", "all_atoms_connected",
    "bond_lengths", "bond_angles", "internal_steric_clash",
    "aromatic_ring_flatness", "double_bond_flatness",
    "internal_energy", "protein-ligand_maximum_distance",
    "minimum_distance_to_protein"
}

# Tools
OST = "/home/rquiroga/anaconda3/envs/runs_n_poses/bin/ost"
PB_VENV_PYTHON = "/home/rquiroga/Datasets/Posebusters/venv/bin/python"


def convert_pdbt_to_sdf(pdbt_file: str, sdf_file: str) -> bool:
    """Convert PDBT to SDF using obabel."""
    try:
        r = subprocess.run(
            ["obabel-25-07", pdbt_file, "-O", sdf_file, "-h"],
            capture_output=True, text=True, timeout=30
        )
        return r.returncode == 0 and os.path.exists(sdf_file) and os.path.getsize(sdf_file) > 100
    except Exception:
        return False


def run_ost_comparison(receptor_cif: str, docked_sdf: str, gt_sdf: str, output_json: str) -> bool:
    """
    Run ost compare-ligand-structures.
    
    Uses the same receptor for both model and reference since we're only
    comparing ligand poses. The docked ligand is provided via -ml flag.
    """
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
                              annotations_row: dict, vinardo_score: float) -> dict:
    """Extract metrics from ost analysis JSON."""
    try:
        with open(json_file) as f:
            result = json.load(f)
    except Exception:
        return None

    if result.get("status") != "SUCCESS":
        return None

    metrics = {
        "target": system_id,
        "method": "vinardock_2vinardo",
        "seed": 1,
        "sample": 1,
        "ranking_score": vinardo_score,
        "ligand_instance_chain": ligand_chain,
        "ligand_is_proper": annotations_row.get("ligand_is_proper", True),
        "model_ligand_ccd_code": annotations_row.get("ligand_ccd_code", ""),
        "model_ligand_smiles": annotations_row.get("ligand_smiles", ""),
        "ligand_ccd_code": annotations_row.get("ligand_ccd_code", ""),
        "model_ligand_chain_lddt_pli": "A",
        "model_ligand_chain_rmsd": "A",
    }

    # Extract lddt_pli
    if "lddt_pli" in result and "assigned_scores" in result["lddt_pli"]:
        for item in result["lddt_pli"]["assigned_scores"]:
            metrics["lddt_pli"] = item.get("score")
            break

    # Extract rmsd
    if "rmsd" in result and "assigned_scores" in result["rmsd"]:
        for item in result["rmsd"]["assigned_scores"]:
            metrics["rmsd"] = item.get("score")
            metrics["lddt_lp"] = item.get("lddt_lp")
            metrics["bb_rmsd"] = item.get("bb_rmsd")
            break

    # iPTM scores - not applicable for docking
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
    """Run PoseBusters on a single docking result.
    
    Converts receptor CIF to PDB since PoseBusters can't read CIF.
    Uses config='redock' to match the AF3 test setup.
    """
    try:
        # Convert receptor CIF to PDB for PoseBusters
        receptor_pdb = receptor_cif.rsplit(".", 1)[0] + ".pdb"
        subprocess.run(
            ["obabel-25-07", receptor_cif, "-O", receptor_pdb],
            capture_output=True, text=True, timeout=30
        )
        
        cmd = [
            PB_VENV_PYTHON, "-c",
            f"""
import json
from posebusters import PoseBusters
bust = PoseBusters(config="redock")
results = bust.bust(
    mol_pred="{model_sdf}",
    mol_true="{gt_sdf}",
    mol_cond="{receptor_pdb}",
)
# Compute pb_success from check columns (matching the paper's criteria)
checks = [
    'sanitization', 'inchi_convertible', 'all_atoms_connected',
    'bond_lengths', 'bond_angles', 'internal_steric_clash',
    'aromatic_ring_flatness', 'double_bond_flatness',
    'internal_energy', 'protein-ligand_maximum_distance',
    'minimum_distance_to_protein'
]
available = [c for c in checks if c in results.columns]
results['pb_success'] = results[available].all(axis=1)
print(results.to_json(orient="records"))
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


def analyze_system(system_id: str, annotations: pd.DataFrame, run_posebusters_flag: bool = True,
                   run_ost_flag: bool = True):
    """Analyze all ligands for a system.

    Parameters:
    - run_posebusters_flag: whether to run PoseBusters for each ligand
    - run_ost_flag: whether to run OST and extract metrics for each ligand
    If run_ost_flag is False, the function will still attempt to run PoseBusters
    (if requested) but will not produce prediction metric rows.
    """
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

        # Check if docking succeeded
        log_file = lig_dir / "log.csv"
        if not log_file.exists() or log_file.stat().st_size < 50:
            continue

        # Extract vinardo score
        vinardo_score = None
        try:
            log_lines = log_file.read_text().strip().split("\n")
            if len(log_lines) >= 2:
                parts_log = log_lines[1].split(",")
                if len(parts_log) >= 3:
                    vinardo_score = float(parts_log[2])
        except Exception:
            pass

        # Find output PDBT
        pdbt_files = [f for f in lig_dir.glob("*.pdbt") if f.stat().st_size > 100]
        if not pdbt_files:
            continue
        vinardo_pdbt = pdbt_files[0]

        # Convert PDBT to SDF for ost -ml
        docked_sdf = lig_dir / f"{chain}_dock.sdf"
        if not docked_sdf.exists() or docked_sdf.stat().st_size < 100:
            if not convert_pdbt_to_sdf(str(vinardo_pdbt), str(docked_sdf)):
                continue

        # Ground truth ligand
        gt_sdf = GT_DIR / system_id / "ligand_files" / f"{chain}.sdf"
        if not gt_sdf.exists():
            continue

        # Run OST comparison and/or PoseBusters depending on flags
        ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
        json_file = ANALYSIS_DIR / f"{system_id}_{chain}.json"

        # Convert and ensure docked_sdf exists (needed for PoseBusters even if OST skipped)
        # (docked_sdf already created above)

        # If OST is requested, run it (unless json exists)
        metrics = None
        if run_ost_flag:
            if not json_file.exists():
                if not run_ost_comparison(str(receptor_cif), str(docked_sdf), str(gt_sdf), str(json_file)):
                    # if OST fails and we're not running PoseBusters, skip this ligand
                    if not run_posebusters_flag:
                        continue
            # Extract metrics if json exists
            if json_file.exists():
                annot_row = system_annot[system_annot["ligand_instance_chain"] == chain]
                if not annot_row.empty:
                    metrics = extract_metrics_from_ost(str(json_file), system_id, chain,
                                                       annot_row.iloc[0].to_dict(), vinardo_score)
                else:
                    metrics = None
        # Run PoseBusters if requested
        if run_posebusters_flag:
            # create temporary RDKit-sanitized SDFs with explicit Hs for PB input (non-destructive)
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

            # Attempt to create a fixed predicted SDF by copying GT valence/H using the fixer script.
            # This is non-destructive: if the fixer fails, we'll fall back to the RDKit-sanitized or original SDF.
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
            # cleanup temp files
            try:
                if tmp_pred and os.path.exists(tmp_pred):
                    os.remove(tmp_pred)
                if tmp_gt and os.path.exists(tmp_gt):
                    os.remove(tmp_gt)
                if tmp_fixed_pred and os.path.exists(tmp_fixed_pred):
                    os.remove(tmp_fixed_pred)
            except Exception:
                pass
            # If metrics are present (we extracted OST), augment with pb_success
            if metrics is not None:
                metrics["pb_success"] = 1.0 if pb_results.get("pb_success", False) else 0.0
            # Build a row for the PoseBusters CSV matching af3.csv-like format
            pb_row = dict(pb_results) if isinstance(pb_results, dict) else {}
            # Ensure a normalized numeric pb_success column exists for downstream plotting
            try:
                pb_row["pb_success"] = 1.0 if pb_results.get("pb_success", False) else 0.0
            except Exception:
                pb_row["pb_success"] = 0.0
            pb_row["system_id"] = system_id
            pb_row.setdefault("seed", (metrics.get("seed", 1) if metrics is not None else 1))
            pb_row.setdefault("sample", (metrics.get("sample", 1) if metrics is not None else 1))
            pb_row.setdefault("ligand_chain", chain)
            pb_row.setdefault("method", "vinardock_2vinardo")
            pb_rows.append(pb_row)
        else:
            if metrics is not None:
                metrics["pb_success"] = -1.0

        if metrics is not None:
            results.append(metrics)

    return results, pb_rows


def make_rdkit_sanitized_sdf(src_path: str) -> str:
    """Create a temporary SDF sanitized with RDKit (explicit Hs).

    If RDKit is not available, fall back to an OpenBabel round-trip with -h.
    Returns path to the temporary SDF file (caller is responsible for cleanup).
    """
    # Try RDKit first
    try:
        from rdkit import Chem
    except Exception:
        Chem = None

    # Fallback to OpenBabel if RDKit isn't available
    if Chem is None:
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
        tmp.close()
        try:
            r = subprocess.run(
                ["obabel-25-07", src_path, "-O", tmp.name, "-h"],
                capture_output=True, text=True, timeout=30
            )
            if r.returncode == 0 and os.path.exists(tmp.name) and os.path.getsize(tmp.name) > 0:
                return tmp.name
        except Exception:
            pass
        if os.path.exists(tmp.name):
            os.remove(tmp.name)
        raise RuntimeError("Neither RDKit available nor OpenBabel round-trip succeeded")

    # Use RDKit to read, sanitize, add explicit Hs and write a new temporary SDF
    suppl = Chem.SDMolSupplier(src_path, removeHs=False)
    mols = [m for m in suppl if m is not None]
    if not mols:
        # try reading removing Hs then add them
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
    # Basic check
    if not os.path.exists(tmp.name) or os.path.getsize(tmp.name) == 0:
        if os.path.exists(tmp.name):
            os.remove(tmp.name)
        raise RuntimeError("Failed to write RDKit-sanitized SDF")
    return tmp.name


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--system-list", type=str,
                        default="/home/rquiroga/github/runs-n-poses/scripts/single_ligand_systems.txt")
    parser.add_argument("--start-index", type=int, default=0)
    parser.add_argument("--max-systems", type=int, default=None)
    parser.add_argument("--skip-posebusters", action="store_true",
                        help="Skip PoseBusters checks")
    parser.add_argument("--resume", action="store_true",
                        help="Resume from existing output CSV and skip already processed systems")
    args = parser.parse_args()

    with open(args.system_list) as f:
        system_ids = [l.strip() for l in f if l.strip()]
    system_ids = system_ids[args.start_index:]
    if args.max_systems:
        system_ids = system_ids[:args.max_systems]

    annotations = pd.read_csv(ANNOTATIONS)
    run_pb = not args.skip_posebusters
    resume = args.resume

    print(f"Analyzing {len(system_ids)} systems (posebusters={'ON' if run_pb else 'OFF'})...")

    # If resuming, read existing Predictions and PoseBusters CSVs to determine what to skip
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

    total_metrics = 0
    # Ensure posebusters_results directory exists
    POSEBUSTERS_CSV.parent.mkdir(parents=True, exist_ok=True)

    for sid in tqdm(system_ids, desc="Analyzing"):
        # Determine per-system resume state
        pred_done = resume and (sid in processed_pred_targets)
        pb_done = resume and (sid in processed_pb_targets)

        # If both parts already done, skip
        if pred_done and pb_done:
            tqdm.write(f"Skipping {sid} (already in {OUTPUT_CSV} and {POSEBUSTERS_CSV})")
            continue

        # Decide which parts to run
        run_ost_flag = not pred_done
        run_pb_for_system = run_pb and (not pb_done)

        metrics, pb_rows = analyze_system(sid, annotations,
                                         run_posebusters_flag=run_pb_for_system,
                                         run_ost_flag=run_ost_flag)

        # If no metrics produced and no PB rows, nothing to write
        if (not metrics) and (not pb_rows):
            continue

        # Append per-system metrics to predictions CSV (only if OST/metrics were run)
        if metrics and run_ost_flag:
            df_sys = pd.DataFrame(metrics)
            write_header = not OUTPUT_CSV.exists()
            try:
                df_sys.to_csv(OUTPUT_CSV, index=False, header=write_header, mode='a')
            except Exception:
                df_sys.to_csv(OUTPUT_CSV, index=False)
            processed_pred_targets.add(sid)

        # Append PoseBusters rows if any; normalize to AF3 canonical columns
        if pb_rows and run_pb_for_system:
            df_pb = pd.DataFrame(pb_rows)
            # Ensure essential metadata columns exist
            for c, default in [("system_id", sid), ("seed", 1), ("sample", 1), ("ligand_chain", None), ("method", "vinardock_2vinardo")]:
                if c not in df_pb.columns:
                    df_pb[c] = default

            # Reindex to canonical AF3 columns; missing columns will be created
            for col in CANONICAL_PB_COLUMNS:
                if col not in df_pb.columns:
                    df_pb[col] = pd.NA
            df_pb = df_pb[CANONICAL_PB_COLUMNS]

            # Fill missing boolean check columns with False (fail) where appropriate
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
        if "lddt_pli" in df_final.columns:
            valid = df_final['lddt_pli'].dropna()
            if len(valid) > 0:
                print(f"  Mean LDDT-PLI: {valid.mean():.3f}")
        if "rmsd" in df_final.columns:
            valid = df_final['rmsd'].dropna()
            if len(valid) > 0:
                print(f"  Mean RMSD: {valid.mean():.3f}")
        if "pb_success" in df_final.columns:
            pb_ok = (df_final['pb_success'] == 1).sum()
            print(f"  PoseBusters pass: {pb_ok}/{len(df_final)} ({100*pb_ok/len(df_final):.1f}%)")
    else:
        print("No new metrics extracted.")


if __name__ == "__main__":
    main()
