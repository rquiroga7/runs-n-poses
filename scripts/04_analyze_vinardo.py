#!/usr/bin/env python
"""
Step 4: Analyze vinardo results and output in the same format as the predictions folder.

This script takes the 2vinardo-mar5 output files and:
1. Scores the predicted poses against ground truth using ost compare-ligand-structures
2. Extracts accuracy metrics (lddt_pli, rmsd, lddt_lp, bb_rmsd, pred_pocket_f1)
3. Outputs a CSV file in the same format as the other prediction methods

The output will be named vinardock_2vinardo.csv and placed in the predictions directory.

Usage:
    python 04_analyze_vinardo.py --vinardo-outputs /path/to/vinardo_outputs \
                                 --ground-truth /path/to/ground_truth \
                                 --annotations /path/to/annotations.csv \
                                 --output /path/to/output/vinardock_2vinardo.csv
"""

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm


def run_ost_comparison(vinardo_output: str, ground_truth_sdf: str, 
                       receptor_cif: str, output_json: str) -> bool:
    """
    Run ost compare-ligand-structures to score the vinardo output.
    
    This uses the same tool as the other methods (AF3, Boltz, etc.) to ensure
    comparable metrics.
    """
    cmd = [
        "ost",
        "compare-ligand-structures",
        "-m", vinardo_output,
        "-rl", ground_truth_sdf,
        "-r", receptor_cif,
        "-o", output_json,
        "--lddt-pli", "--rmsd", "--lddt-pli-amc"
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120
        )
        
        if result.returncode != 0:
            print(f"Warning: ost compare failed for {vinardo_output}:")
            print(f"  stderr: {result.stderr[:300]}")
            return False
        
        return os.path.exists(output_json)
        
    except subprocess.TimeoutExpired:
        print(f"Warning: ost compare timed out for {vinardo_output}")
        return False
    except Exception as e:
        print(f"Error running ost compare for {vinardo_output}: {e}")
        return False


def extract_metrics_from_json(json_file: str, system_id: str, ligand_chain: str,
                               annotations_row: dict) -> dict:
    """
    Extract metrics from the ost compare-ligand-structures JSON output.
    
    Returns a dictionary with the same columns as the prediction CSV files.
    """
    try:
        with open(json_file) as f:
            result = json.load(f)
    except Exception as e:
        print(f"Error reading JSON {json_file}: {e}")
        return None
    
    # Initialize result dictionary with required columns
    metrics = {
        "target": system_id,
        "method": "vinardock_2vinardo",
        "seed": 1,  # Vinardo doesn't use seeds
        "sample": 1,  # Single sample
        "ranking_score": None,  # Vinardo scoring
        "ligand_instance_chain": ligand_chain,
        "ligand_is_proper": annotations_row.get("ligand_is_proper", True),
    }
    
    # Extract lddt_pli metrics
    if "lddt_pli" in result and "assigned_scores" in result["lddt_pli"]:
        for item in result["lddt_pli"]["assigned_scores"]:
            if item.get("reference_ligand") and ligand_chain in item.get("reference_ligand", ""):
                metrics["lddt_pli"] = item.get("score")
                metrics["model_ligand_chain_lddt_pli"] = item.get("model_ligand", "").split(".")[0]
                break
    
    # Extract RMSD metrics
    if "rmsd" in result and "assigned_scores" in result["rmsd"]:
        for item in result["rmsd"]["assigned_scores"]:
            if item.get("reference_ligand") and ligand_chain in item.get("reference_ligand", ""):
                metrics["rmsd"] = item.get("score")
                metrics["lddt_lp"] = item.get("lddt_lp")
                metrics["bb_rmsd"] = item.get("bb_rmsd")
                metrics["model_ligand_chain_rmsd"] = item.get("model_ligand", "").split(".")[0]
                break
    
    # Add ligand information from annotations
    metrics["model_ligand_ccd_code"] = annotations_row.get("ligand_ccd_code", "")
    metrics["model_ligand_smiles"] = annotations_row.get("ligand_smiles", "")
    metrics["ligand_ccd_code"] = annotations_row.get("ligand_ccd_code", "")
    
    # iPTM scores - not applicable for vinardo, set to None
    for col in [
        "prot_lig_chain_iptm_average_lddt_pli",
        "prot_lig_chain_iptm_min_lddt_pli",
        "prot_lig_chain_iptm_max_lddt_pli",
        "lig_prot_chain_iptm_average_lddt_pli",
        "lig_prot_chain_iptm_min_lddt_pli",
        "lig_prot_chain_iptm_max_lddt_pli",
        "prot_lig_chain_iptm_average_rmsd",
        "prot_lig_chain_iptm_min_rmsd",
        "prot_lig_chain_iptm_max_rmsd",
        "lig_prot_chain_iptm_average_rmsd",
        "lig_prot_chain_iptm_min_rmsd",
        "lig_prot_chain_iptm_max_rmsd",
    ]:
        metrics[col] = None
    
    # Pocket prediction metrics - not applicable for vinardo
    metrics["pred_pocket_tp"] = None
    metrics["pred_pocket_fp"] = None
    metrics["pred_pocket_fn"] = None
    metrics["pred_pocket_precision"] = None
    metrics["pred_pocket_recall"] = None
    metrics["pred_pocket_f1"] = None
    
    return metrics


def extract_vinardo_score(output_file: str) -> float:
    """
    Extract the ranking score from vinardo output.
    
    The vinardo output file contains the docking score.
    """
    try:
        with open(output_file) as f:
            content = f.read()
        
        # Look for REMARK lines with scores
        for line in content.split("\n"):
            if line.startswith("REMARK VINA RESULT"):
                # Format: REMARK VINA RESULT:   -9.2     0.000     0.000
                parts = line.split(":")
                if len(parts) > 1:
                    scores = parts[1].strip().split()
                    if scores:
                        return float(scores[0])
    except Exception:
        pass
    
    return None


def process_system_output(system_id: str, ligand_chain: str, 
                          vinardo_output_dir: Path, ground_truth_dir: Path,
                          annotations: pd.DataFrame, analysis_output_dir: Path) -> dict:
    """
    Process a single vinardo output and extract metrics.
    """
    # Get vinardo output file
    vinardo_output_file = vinardo_output_dir / system_id / f"{system_id}_{ligand_chain}_output.pdbqt"
    
    if not vinardo_output_file.exists():
        print(f"Warning: Vinardo output not found for {system_id} {ligand_chain}")
        return None
    
    # Get ground truth files
    system_gt_dir = ground_truth_dir / system_id
    ground_truth_sdf = system_gt_dir / "ligand_files" / f"{ligand_chain}.sdf"
    receptor_cif = system_gt_dir / "receptor.cif"
    
    if not ground_truth_sdf.exists():
        print(f"Warning: Ground truth SDF not found for {system_id} {ligand_chain}")
        return None
    
    if not receptor_cif.exists():
        print(f"Warning: Receptor CIF not found for {system_id}")
        return None
    
    # Get annotation row for this system
    system_annotations = annotations[
        (annotations["system_id"] == system_id) & 
        (annotations["ligand_instance_chain"] == ligand_chain)
    ]
    
    if system_annotations.empty:
        print(f"Warning: No annotations for {system_id} {ligand_chain}")
        return None
    
    annotations_row = system_annotations.iloc[0].to_dict()
    
    # Run ost comparison
    analysis_output_dir.mkdir(parents=True, exist_ok=True)
    output_json = analysis_output_dir / f"{system_id}_{ligand_chain}.json"
    
    if not run_ost_comparison(
        str(vinardo_output_file),
        str(ground_truth_sdf),
        str(receptor_cif),
        str(output_json)
    ):
        return None
    
    # Extract metrics
    metrics = extract_metrics_from_json(
        str(output_json),
        system_id,
        ligand_chain,
        annotations_row
    )
    
    # Add vinardo scoring as ranking_score
    if metrics:
        metrics["ranking_score"] = extract_vinardo_score(str(vinardo_output_file))
    
    return metrics


def main():
    parser = argparse.ArgumentParser(description="Analyze vinardo results and output in predictions format")
    parser.add_argument(
        "--vinardo-outputs",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_outputs",
        help="Path to vinardo outputs directory"
    )
    parser.add_argument(
        "--ground-truth",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth",
        help="Path to ground_truth directory"
    )
    parser.add_argument(
        "--annotations",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/annotations.csv",
        help="Path to annotations.csv file"
    )
    parser.add_argument(
        "--inputs-json",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/inputs.json",
        help="Path to inputs.json file"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/predictions/vinardock_2vinardo.csv",
        help="Path to output CSV file"
    )
    parser.add_argument(
        "--analysis-dir",
        type=str,
        default="/home/rquiroga/github/runs-n-poses/examples/analysis/vinardock_2vinardo",
        help="Path to store intermediate analysis JSON files (follows examples/analysis/ pattern)"
    )
    parser.add_argument(
        "--system-id",
        type=str,
        default=None,
        help="Process only a specific system ID (optional)"
    )
    
    args = parser.parse_args()

    # Load annotations
    if not os.path.exists(args.annotations):
        print(f"Error: Annotations file not found: {args.annotations}")
        sys.exit(1)

    annotations = pd.read_csv(args.annotations)

    # Load inputs.json (following the pattern from extract_scores.ipynb)
    if not os.path.exists(args.inputs_json):
        print(f"Warning: inputs.json not found at {args.inputs_json}")
        input_data = {}
    else:
        with open(args.inputs_json, 'r') as f:
            input_data = json.load(f)
        print(f"Loaded inputs.json with {len(input_data)} systems")

    vinardo_outputs_dir = Path(args.vinardo_outputs)
    ground_truth_dir = Path(args.ground_truth)
    analysis_output_dir = Path(args.analysis_dir)
    analysis_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get list of systems to process
    if args.system_id:
        system_ids = [args.system_id]
    else:
        # Get unique system-ligand pairs from annotations
        system_ligand_pairs = annotations[["system_id", "ligand_instance_chain"]].drop_duplicates()
        system_ids = sorted(system_ligand_pairs["system_id"].unique())
    
    print(f"Processing {len(system_ids)} systems...")
    
    all_metrics = []
    success_count = 0
    fail_count = 0
    
    for system_id in tqdm(system_ids, desc="Analyzing vinardo outputs"):
        # Get all ligand chains for this system
        system_ligands = annotations[
            annotations["system_id"] == system_id
        ]["ligand_instance_chain"].unique()
        
        for ligand_chain in system_ligands:
            metrics = process_system_output(
                system_id,
                ligand_chain,
                vinardo_outputs_dir,
                ground_truth_dir,
                annotations,
                analysis_output_dir
            )
            
            if metrics:
                all_metrics.append(metrics)
                success_count += 1
            else:
                fail_count += 1
    
    # Create output DataFrame
    if all_metrics:
        results_df = pd.DataFrame(all_metrics)
        
        # Ensure output directory exists
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Save to CSV
        results_df.to_csv(output_path, index=False)
        
        print(f"\nDone! Success: {success_count}, Failed: {fail_count}")
        print(f"Results saved to: {output_path}")
        print(f"\nSummary:")
        print(f"  Total systems: {results_df['target'].nunique()}")
        print(f"  Total predictions: {len(results_df)}")
        if "lddt_pli" in results_df.columns:
            print(f"  Mean LDDT-PLI: {results_df['lddt_pli'].mean():.3f}")
        if "rmsd" in results_df.columns:
            print(f"  Mean RMSD: {results_df['rmsd'].mean():.3f}")
    else:
        print(f"\nWarning: No metrics extracted. Output file not created.")
        print(f"Success: {success_count}, Failed: {fail_count}")


if __name__ == "__main__":
    main()
