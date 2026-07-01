#!/usr/bin/env python
"""
Step 3c: Prepare ligand PDBT and PDBQT files for all systems.

Reads the proper ligand SDF from ground_truth/<sys_id>/ligand_files/<chain>.sdf
and converts to both PDBT and PDBQT formats for docking.

Usage:
    python prepare_symmetry_ligands.py \
        --ground-truth-dir runs-n-poses-datasets/ground_truth \
        --output-pdbt runs-n-poses-datasets/symmetry_ligands_pdbt \
        --output-pdbqt runs-n-poses-datasets/symmetry_ligands_pdbqt \
        --system-info scripts/systems_for_symmetry_docking.csv
"""

import argparse
import csv
import os
import subprocess
from pathlib import Path

OBABEL = "obabel-25-07"


def convert_sdf(sdf_path, output_path):
    """Convert SDF to PDBT/PDBQT using obabel (format inferred from extension)."""
    cmd = [OBABEL, str(sdf_path), "-O", str(output_path), "-p", "7.4"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            print(f"    obabel stderr: {result.stderr[:200]}")
            return False
        return os.path.exists(output_path) and os.path.getsize(output_path) > 0
    except Exception as e:
        print(f"    obabel exception: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(description="Prepare ligand PDBT and PDBQT files")
    parser.add_argument("--ground-truth-dir", required=True)
    parser.add_argument("--output-pdbt", required=True)
    parser.add_argument("--output-pdbqt", required=True)
    parser.add_argument("--system-info", required=True)
    parser.add_argument("--system-id", type=str, default=None)
    args = parser.parse_args()

    gt_dir = Path(args.ground_truth_dir)
    out_pdbt = Path(args.output_pdbt)
    out_pdbqt = Path(args.output_pdbqt)

    with open(args.system_info) as f:
        reader = csv.DictReader(f)
        systems = [(r["system_id"], r["proper_ligand_chain"]) for r in reader]
    if args.system_id:
        systems = [s for s in systems if s[0] == args.system_id]

    ok_pdbt, ok_pdbqt, fail = 0, 0, 0
    for sys_id, chain in systems:
        sdf_path = gt_dir / sys_id / "ligand_files" / f"{chain}.sdf"
        if not sdf_path.exists():
            print(f"  WARN: {sdf_path} not found")
            fail += 1
            continue

        # PDBT
        pdbt_dir = out_pdbt / sys_id
        pdbt_dir.mkdir(parents=True, exist_ok=True)
        pdbt_path = pdbt_dir / f"{chain}.pdbt"
        if not (pdbt_path.exists() and pdbt_path.stat().st_size > 0):
            if convert_sdf(str(sdf_path), str(pdbt_path)):
                ok_pdbt += 1
            else:
                print(f"  FAIL pdbt: {sys_id} {chain}")
                fail += 1
        else:
            ok_pdbt += 1

        # PDBQT
        pdbqt_dir = out_pdbqt / sys_id
        pdbqt_dir.mkdir(parents=True, exist_ok=True)
        pdbqt_path = pdbqt_dir / f"{chain}.pdbqt"
        if not (pdbqt_path.exists() and pdbqt_path.stat().st_size > 0):
            if convert_sdf(str(sdf_path), str(pdbqt_path)):
                ok_pdbqt += 1
            else:
                print(f"  FAIL pdbqt: {sys_id} {chain}")
                fail += 1
        else:
            ok_pdbqt += 1

    print(f"Ligands: PDBT={ok_pdbt}, PDBQT={ok_pdbqt}, Fail={fail}")


if __name__ == "__main__":
    main()
