#!/usr/bin/env python
"""
Step 3a: Convert symmetry-corrected receptor PDB files to PDBT format for Vinardo.

Reads symmetry_corrected/<sys_id>/<sys_id>_receptor_symm.pdb and writes
symmetry_receptors_pdbt/<sys_id>/<sys_id>_receptor.pdbt via obabel.

Usage:
    python prepare_receptor_pdbt.py \
        --input-dir runs-n-poses-datasets/symmetry_corrected \
        --output-dir runs-n-poses-datasets/symmetry_receptors_pdbt \
        --system-list scripts/systems_for_symmetry_docking.txt
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

OBABEL = "obabel-25-07"


def convert_pdb_to_pdbt(pdb_path, output_path):
    """Convert PDB to PDBT using obabel with -h -xc -xr flags."""
    cmd = [
        OBABEL,
        pdb_path,
        "-O", output_path,
        "-p 7.4",
        "-xc",
        "-xr",
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            return False, result.stderr
        if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
            return True, None
        return False, "Output file empty"
    except subprocess.TimeoutExpired:
        return False, "Timeout"
    except Exception as e:
        return False, str(e)


def main():
    parser = argparse.ArgumentParser(description="Convert symm PDB to PDBT for Vinardo")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--system-list", required=True)
    parser.add_argument("--system-id", type=str, default=None)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    with open(args.system_list) as f:
        system_ids = [l.strip() for l in f if l.strip()]
    if args.system_id:
        system_ids = [s for s in system_ids if s == args.system_id]

    ok, fail, skip = 0, 0, 0
    for sys_id in system_ids:
        pdb_path = input_dir / sys_id / f"{sys_id}_receptor_symm.pdb"
        out_dir = output_dir / sys_id
        out_path = out_dir / f"{sys_id}_receptor.pdbt"
        out_dir.mkdir(parents=True, exist_ok=True)

        if out_path.exists() and out_path.stat().st_size > 0:
            skip += 1
            continue

        if not pdb_path.exists():
            print(f"  WARN: {pdb_path} not found")
            fail += 1
            continue

        success, err = convert_pdb_to_pdbt(str(pdb_path), str(out_path))
        if success:
            ok += 1
        else:
            print(f"  FAIL: {sys_id}: {err}")
            fail += 1

    print(f"Receptor PDBT: OK={ok}, Skip={skip}, Fail={fail}")


if __name__ == "__main__":
    main()
