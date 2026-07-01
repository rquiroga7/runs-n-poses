#!/usr/bin/env python
"""
Step 3b: Convert symmetry-corrected receptor PDB files to PDBQT format for Vina.

Reads symmetry_corrected/<sys_id>/<sys_id>_receptor_symm.pdb and writes
symmetry_receptors_pdbqt/<sys_id>/<sys_id>_receptor.pdbqt via obabel.

If the obabel output contains ROOT/ENDROOT tags (which Vina rejects), a
round-trip PDBQT -> PDB -> PDBQT cleanup is performed automatically.

Usage:
    python prepare_receptor_pdbqt.py \
        --input-dir runs-n-poses-datasets/symmetry_corrected \
        --output-dir runs-n-poses-datasets/symmetry_receptors_pdbqt \
        --system-list scripts/systems_for_symmetry_docking.txt
"""

import argparse
import os
import subprocess
import tempfile
from pathlib import Path

OBABEL = "obabel-25-07"


def has_root_tags(pdbqt_path):
    """Check if a PDBQT file contains ROOT/ENDROOT tags."""
    try:
        with open(pdbqt_path) as f:
            content = f.read()
        return "ROOT" in content or "ENDROOT" in content
    except Exception:
        return False


def clean_root_tags(pdbqt_path):
    """Round-trip PDBQT -> PDB -> PDBQT to remove ROOT tags.
    
    Operates in-place on the file.
    """
    stem = str(pdbqt_path).replace(".pdbqt", "")
    try:
        with tempfile.TemporaryDirectory() as tmp:
            roundtrip_pdb = os.path.join(tmp, "roundtrip.pdb")
            roundtrip_pdbqt = os.path.join(tmp, "roundtrip.pdbqt")

            r1 = subprocess.run(
                [OBABEL, str(pdbqt_path), "-O", roundtrip_pdb],
                capture_output=True, text=True, timeout=60
            )
            if r1.returncode != 0 or not os.path.exists(roundtrip_pdb):
                return False

            r2 = subprocess.run(
                [OBABEL, roundtrip_pdb, "-O", roundtrip_pdbqt, "-h", "-xc", "-xr"],
                capture_output=True, text=True, timeout=60
            )
            if r2.returncode != 0 or not os.path.exists(roundtrip_pdbqt):
                return False

            import shutil
            shutil.copy2(roundtrip_pdbqt, str(pdbqt_path))
            return True
    except Exception:
        return False


def convert_pdb_to_pdbqt(pdb_path, output_path):
    """Convert PDB to PDBQT using obabel. Returns success flag."""
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
        if not (os.path.exists(output_path) and os.path.getsize(output_path) > 0):
            return False, "Output file empty"
        # Clean ROOT tags if present
        if has_root_tags(output_path):
            if not clean_root_tags(output_path):
                return False, "ROOT cleanup failed"
            if has_root_tags(output_path):
                return False, "ROOT tags still present after cleanup"
        return True, None
    except subprocess.TimeoutExpired:
        return False, "Timeout"
    except Exception as e:
        return False, str(e)


def main():
    parser = argparse.ArgumentParser(description="Convert symm PDB to PDBQT for Vina")
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
        out_path = out_dir / f"{sys_id}_receptor.pdbqt"
        out_dir.mkdir(parents=True, exist_ok=True)

        if out_path.exists() and out_path.stat().st_size > 0:
            skip += 1
            continue

        if not pdb_path.exists():
            print(f"  WARN: {pdb_path} not found")
            fail += 1
            continue

        success, err = convert_pdb_to_pdbqt(str(pdb_path), str(out_path))
        if success:
            ok += 1
        else:
            print(f"  FAIL: {sys_id}: {err}")
            fail += 1

    print(f"Receptor PDBQT: OK={ok}, Skip={skip}, Fail={fail}")


if __name__ == "__main__":
    main()
