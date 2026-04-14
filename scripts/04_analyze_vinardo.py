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

GT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth")
OUTPUT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_outputs")
ANALYSIS_DIR = Path("/home/rquiroga/github/runs-n-poses/examples/analysis/vinardock_2vinardo")
ANNOTATIONS = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/annotations.csv")
OUTPUT_CSV = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/predictions/vinardock_2vinardo.csv")
POSEBUSTERS_CSV = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/posebusters_results/vinardock_2vinardo.csv")

# Tools
OST = "/home/rquiroga/anaconda3/envs/runs_n_poses/bin/ost"
PB_VENV_PYTHON = "/home/rquiroga/Datasets/Posebusters/venv/bin/python"


def convert_pdbt_to_sdf(pdbt_file: str, sdf_file: str) -> bool:
    """Convert PDBT to SDF using obabel."""
    try:
        r = subprocess.run(
            ["obabel-25-07", pdbt_file, "-O", sdf_file],
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


def analyze_system(system_id: str, annotations: pd.DataFrame, run_posebusters_flag: bool = True) -> list:
    """Analyze all ligands for a system."""
    system_annot = annotations[annotations["system_id"] == system_id]
    if system_annot.empty:
        return []

    results = []
    system_out = OUTPUT_DIR / system_id
    if not system_out.exists():
        return []

    receptor_cif = GT_DIR / system_id / "receptor.cif"
    if not receptor_cif.exists():
        return []

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

        # Run ost comparison (ligand-only comparison against ground truth)
        ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
        json_file = ANALYSIS_DIR / f"{system_id}_{chain}.json"

        if not json_file.exists():
            if not run_ost_comparison(str(receptor_cif), str(docked_sdf), str(gt_sdf), str(json_file)):
                continue

        # Extract metrics
        annot_row = system_annot[system_annot["ligand_instance_chain"] == chain]
        if annot_row.empty:
            continue

        metrics = extract_metrics_from_ost(str(json_file), system_id, chain,
                                           annot_row.iloc[0].to_dict(), vinardo_score)
        if metrics:
            # Run PoseBusters
            if run_posebusters_flag:
                pb_results = run_posebusters(str(docked_sdf), str(gt_sdf), str(receptor_cif))
                metrics["pb_success"] = 1.0 if pb_results.get("pb_success", False) else 0.0
            else:
                metrics["pb_success"] = -1.0
            results.append(metrics)

    return results


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

    # If resuming, read existing output CSV to find already processed systems
    processed_targets = set()
    if resume and OUTPUT_CSV.exists():
        try:
            existing = pd.read_csv(OUTPUT_CSV)
            if 'target' in existing.columns:
                processed_targets = set(existing['target'].astype(str).unique())
        except Exception:
            processed_targets = set()

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    total_metrics = 0
    for sid in tqdm(system_ids, desc="Analyzing"):
        if resume and sid in processed_targets:
            tqdm.write(f"Skipping {sid} (already present in {OUTPUT_CSV})")
            continue

        metrics = analyze_system(sid, annotations, run_posebusters_flag=run_pb)
        if not metrics:
            continue

        # Append per-system results to CSV to save progress incrementally
        df_sys = pd.DataFrame(metrics)
        write_header = not OUTPUT_CSV.exists()
        try:
            df_sys.to_csv(OUTPUT_CSV, index=False, header=write_header, mode='a')
        except Exception:
            # fallback: try writing normally
            df_sys.to_csv(OUTPUT_CSV, index=False)

        total_metrics += len(df_sys)
        processed_targets.add(sid)

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
