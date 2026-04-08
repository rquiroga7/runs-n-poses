#!/usr/bin/env python
"""
Step 2: Prepare ligand PDBQT files for 2vinardo-mar5.

This script extracts ligand SMILES from annotations.csv for the specified 
ligand_instance_chain and converts them to PDBQT format using obabel-25-07.

Usage:
    python 02_prepare_ligand_pdbqt.py --annotations /path/to/annotations.csv \
                                      --ground-truth-dir /path/to/ground_truth \
                                      --output-dir /path/to/output
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm


def convert_smiles_to_pdbqt(smiles: str, output_file: str, obabel_path: str = "obabel-25-07") -> bool:
    """
    Convert a SMILES string to PDBQT format using obabel.

    -d: add hydrogens (with correct protonation state)
    Note: We do NOT use -xc or -xr for ligands (unlike receptors).
    """
    cmd = [
        obabel_path,
        f"-:{smiles}",
        "-osmi",  # input is SMILES
        "-O", output_file,
        "--gen3d",  # Generate 3D coordinates
        "-d"  # add hydrogens
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120
        )
        if result.returncode != 0:
            print(f"Warning: obabel failed for SMILES {smiles[:50]}...: {result.stderr}")
            return False
        return os.path.exists(output_file) and os.path.getsize(output_file) > 0
    except subprocess.TimeoutExpired:
        print(f"Warning: obabel timed out for SMILES {smiles[:50]}...")
        return False
    except Exception as e:
        print(f"Error converting SMILES {smiles[:50]}...: {e}")
        return False


def extract_ligand_from_sdf(sdf_file: str, output_file: str, obabel_path: str = "obabel-25-07") -> bool:
    """
    Extract ligand from SDF file and convert to PDBQT.
    
    This is an alternative approach if we want to use the ground truth ligand pose.
    """
    cmd = [
        obabel_path,
        sdf_file,
        "-O", output_file
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120
        )
        if result.returncode != 0:
            print(f"Warning: obabel failed for {sdf_file}: {result.stderr}")
            return False
        return os.path.exists(output_file) and os.path.getsize(output_file) > 0
    except subprocess.TimeoutExpired:
        print(f"Warning: obabel timed out for {sdf_file}")
        return False
    except Exception as e:
        print(f"Error converting {sdf_file}: {e}")
        return False


def process_system(system_id: str, annotations: pd.DataFrame, ground_truth_dir: Path, 
                   output_dir: Path, obabel_path: str) -> bool:
    """
    Process a single system and create its ligand PDBQT file(s).
    
    For each ligand_instance_chain in the annotations for this system,
    create a PDBQT file from the SMILES string.
    """
    system_annotations = annotations[annotations["system_id"] == system_id]
    
    if system_annotations.empty:
        print(f"Warning: No annotations found for {system_id}")
        return False
    
    output_system_dir = output_dir / system_id
    output_system_dir.mkdir(parents=True, exist_ok=True)
    
    success = False
    
    for idx, row in system_annotations.iterrows():
        ligand_chain = row["ligand_instance_chain"]
        smiles = row.get("ligand_smiles")
        
        if pd.isna(smiles) or not smiles:
            print(f"Warning: No SMILES for {system_id}, ligand {ligand_chain}")
            # Try to use the SDF file from ground truth instead
            sdf_file = ground_truth_dir / system_id / "ligand_files" / f"{ligand_chain}.sdf"
            if sdf_file.exists():
                output_file = output_system_dir / f"{ligand_chain}.pdbqt"
                if extract_ligand_from_sdf(str(sdf_file), str(output_file), obabel_path):
                    success = True
            continue
        
        output_file = output_system_dir / f"{ligand_chain}.pdbqt"
        if convert_smiles_to_pdbqt(smiles, str(output_file), obabel_path):
            success = True
        else:
            # Fallback: try to use the SDF file
            sdf_file = ground_truth_dir / system_id / "ligand_files" / f"{ligand_chain}.sdf"
            if sdf_file.exists():
                if extract_ligand_from_sdf(str(sdf_file), str(output_file), obabel_path):
                    success = True
    
    return success


def main():
    parser = argparse.ArgumentParser(description="Prepare ligand PDBQT files for 2vinardo-mar5")
    parser.add_argument(
        "--annotations",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/annotations.csv",
        help="Path to annotations.csv file"
    )
    parser.add_argument(
        "--ground-truth-dir",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth",
        help="Path to the ground_truth directory"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_inputs/ligands",
        help="Path to output directory for ligand PDBQT files"
    )
    parser.add_argument(
        "--obabel-path",
        type=str,
        default="obabel-25-07",
        help="Path to obabel-25-07 executable"
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
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    ground_truth_dir = Path(args.ground_truth_dir)
    
    # Get list of systems to process
    if args.system_id:
        system_ids = [args.system_id]
    else:
        system_ids = annotations["system_id"].unique()
        system_ids = sorted(system_ids)
    
    print(f"Processing {len(system_ids)} systems...")
    
    success_count = 0
    fail_count = 0
    
    for system_id in tqdm(system_ids, desc="Preparing ligands"):
        if process_system(system_id, annotations, ground_truth_dir, output_dir, args.obabel_path):
            success_count += 1
        else:
            fail_count += 1
    
    print(f"\nDone! Success: {success_count}, Failed: {fail_count}")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()
