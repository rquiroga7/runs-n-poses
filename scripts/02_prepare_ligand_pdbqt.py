#!/usr/bin/env python
"""
Step 2: Prepare ligand PDBT files for vinardock.

This script takes ligand SDF files from the ground_truth directory and converts 
them to PDBT format using obabel-25-07 with -h flag.

Usage:
    python 02_prepare_ligand_pdbqt.py --ground-truth-dir /path/to/ground_truth \
                                      --output-dir /path/to/output
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

from tqdm import tqdm


def convert_sdf_to_pdbt(sdf_file: str, output_file: str, obabel_path: str = "obabel-25-07") -> bool:
    """
    Convert SDF file to PDBT format using obabel with -h flag.
    
    Uses the ground truth crystal structure coordinates.
    """
    cmd = [
        obabel_path,
        sdf_file,
        "-O", output_file,
        "-p 7.4"
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


def process_system(system_id: str, ground_truth_dir: Path,
                   output_dir: Path, obabel_path: str) -> bool:
    """
    Process a single system and create its ligand PDBT file(s).
    Uses ground truth SDF files directly for proper crystal coordinates.
    """
    system_gt_dir = ground_truth_dir / system_id / "ligand_files"
    
    if not system_gt_dir.exists():
        print(f"Warning: No ligand_files dir for {system_id}")
        return False

    output_system_dir = output_dir / system_id
    output_system_dir.mkdir(parents=True, exist_ok=True)

    success = False

    for sdf_file in sorted(system_gt_dir.glob("*.sdf")):
        ligand_chain = sdf_file.stem
        output_file = output_system_dir / f"{ligand_chain}.pdbt"
        
        if convert_sdf_to_pdbt(str(sdf_file), str(output_file), obabel_path):
            success = True

    return success


def main():
    parser = argparse.ArgumentParser(description="Prepare ligand PDBT files for vinardock")
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
        help="Path to output directory for ligand PDBT files"
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

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ground_truth_dir = Path(args.ground_truth_dir)
    
    # Get system list from ground truth directory
    if args.system_id:
        system_ids = [args.system_id]
    else:
        system_ids = [d.name for d in ground_truth_dir.iterdir() if d.is_dir()]
        system_ids = sorted(system_ids)

    print(f"Processing {len(system_ids)} systems...")

    success_count = 0
    fail_count = 0

    for system_id in tqdm(system_ids, desc="Preparing ligands"):
        if process_system(system_id, ground_truth_dir, output_dir, args.obabel_path):
            success_count += 1
        else:
            fail_count += 1

    print(f"\nDone! Success: {success_count}, Failed: {fail_count}")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()
