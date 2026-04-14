#!/usr/bin/env python
"""
Step 1: Prepare receptor PDBT files for vinardock.

This script takes the receptor.cif files from the ground_truth directory and converts them
to PDBT format using obabel-25-07 with -d -xc -xr flags.

Usage:
    python 01_prepare_receptor_pdbqt.py --ground-truth-dir /path/to/ground_truth \
                                        --output-dir /path/to/output
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm


def find_receptor_and_cofactors(system_dir: Path) -> list:
    """Find receptor.cif in the system directory."""
    receptor_cif = system_dir / "receptor.cif"
    if not receptor_cif.exists():
        return []
    return [str(receptor_cif)]


def convert_to_pdbt(input_file: str, output_file: str, obabel_path: str = "obabel-25-07") -> bool:
    """
    Convert a CIF/PDB file to PDBT format using obabel.

    -h:   add hydrogens
    -xc:  remove charges
    -xr:  remove residues (keeps only organic/coordination compounds)
    """
    cmd = [
        obabel_path,
        input_file,
        "-O", output_file,
        "-h",   # add hydrogens
        "-xc",  # remove charges
        "-xr"   # remove residues
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120
        )
        if result.returncode != 0:
            print(f"Warning: obabel failed for {input_file}: {result.stderr}")
            return False
        return os.path.exists(output_file) and os.path.getsize(output_file) > 0
    except subprocess.TimeoutExpired:
        print(f"Warning: obabel timed out for {input_file}")
        return False
    except Exception as e:
        print(f"Error converting {input_file}: {e}")
        return False


def process_system(system_id: str, ground_truth_dir: Path, output_dir: Path, 
                   obabel_path: str) -> bool:
    """Process a single system and create its receptor PDBT file."""
    system_dir = ground_truth_dir / system_id
    output_system_dir = output_dir / system_id
    output_system_dir.mkdir(parents=True, exist_ok=True)
    
    receptor_files = find_receptor_and_cofactors(system_dir)
    
    if not receptor_files:
        print(f"Warning: No receptor files found for {system_id}")
        return False
    
    success = True
    for receptor_file in receptor_files:
        output_file = output_system_dir / f"{system_id}_receptor.pdbt"
        if not convert_to_pdbt(receptor_file, str(output_file), obabel_path):
            success = False
    
    return success


def main():
    parser = argparse.ArgumentParser(description="Prepare receptor PDBT files for vinardock")
    parser.add_argument(
        "--ground-truth-dir",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth",
        help="Path to the ground_truth directory"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_inputs/receptors",
        help="Path to output directory for PDBT files"
    )
    parser.add_argument(
        "--annotations",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/annotations.csv",
        help="Path to annotations.csv file"
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
    
    ground_truth_dir = Path(args.ground_truth_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if args.system_id:
        system_ids = [args.system_id]
    else:
        system_ids = [d.name for d in ground_truth_dir.iterdir() if d.is_dir()]
        system_ids.sort()
    
    print(f"Processing {len(system_ids)} systems...")
    
    success_count = 0
    fail_count = 0
    
    for system_id in tqdm(system_ids, desc="Preparing receptors"):
        if process_system(system_id, ground_truth_dir, output_dir, args.obabel_path):
            success_count += 1
        else:
            fail_count += 1
    
    print(f"\nDone! Success: {success_count}, Failed: {fail_count}")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()
