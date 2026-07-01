#!/usr/bin/env python
"""
Step 2 (VINA): Prepare ligand PDBQT files for AutoDock Vina.

This script converts ground-truth ligand SDF files to PDBQT using obabel-25-07.

Usage:
    python 02_prepare_ligand_pdbqt_vina.py --ground-truth-dir /path/to/ground_truth \
                                          --output-dir /path/to/output
"""

import argparse
import os
import subprocess
from pathlib import Path
from tqdm import tqdm


def convert_sdf_to_pdbqt(sdf_file: str, output_file: str, obabel_path: str = "obabel-25-07") -> bool:
    cmd = [
        obabel_path,
        sdf_file,
        "-O", output_file,
        "-p 7.4",
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


def process_system(system_id: str, ground_truth_dir: Path, output_dir: Path, obabel_path: str) -> bool:
    system_gt_dir = ground_truth_dir / system_id / "ligand_files"
    if not system_gt_dir.exists():
        print(f"Warning: No ligand_files dir for {system_id}")
        return False

    output_system_dir = output_dir / system_id
    output_system_dir.mkdir(parents=True, exist_ok=True)

    success = False
    for sdf_file in sorted(system_gt_dir.glob("*.sdf")):
        ligand_chain = sdf_file.stem
        output_file = output_system_dir / f"{ligand_chain}.pdbqt"
        if convert_sdf_to_pdbqt(str(sdf_file), str(output_file), obabel_path):
            success = True

    return success


def main():
    parser = argparse.ArgumentParser(description="Prepare ligand PDBQT files for Vina")
    parser.add_argument(
        "--ground-truth-dir",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/vina_inputs/ligands",
    )
    parser.add_argument(
        "--obabel-path",
        type=str,
        default="obabel-25-07",
    )
    parser.add_argument("--system-id", type=str, default=None)

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ground_truth_dir = Path(args.ground_truth_dir)
    if args.system_id:
        system_ids = [args.system_id]
    else:
        system_ids = [d.name for d in ground_truth_dir.iterdir() if d.is_dir()]
        system_ids = sorted(system_ids)

    success_count = 0
    fail_count = 0

    for system_id in tqdm(system_ids, desc="Preparing ligands (Vina)"):
        if process_system(system_id, ground_truth_dir, output_dir, args.obabel_path):
            success_count += 1
        else:
            fail_count += 1

    print(f"\nDone! Success: {success_count}, Failed: {fail_count}")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()
