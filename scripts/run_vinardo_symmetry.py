#!/usr/bin/env python
"""
Step 4b: Run VinardoCK (2vinardo scoring) on symmetry-corrected receptors.

Adapted from 03_run_vinardo.py for the symmetry-corrected docking pipeline.

Usage:
    python run_vinardo_symmetry.py \
        --receptor-dir runs-n-poses-datasets/symmetry_receptors_pdbt \
        --ligand-dir runs-n-poses-datasets/symmetry_ligands_pdbt \
        --output-dir runs-n-poses-datasets/vinardo_outputs_symmetry \
        --system-list scripts/systems_for_symmetry_docking.txt
"""

import argparse
import json
import os
import subprocess
import time
from pathlib import Path

VINARDO = "vinardock-26-04"


def run_vinardock(receptor_file, ligand_file, output_dir, sys_id, lig_chain, threads=8):
    """Run Vinardo docking for a single receptor-ligand pair."""
    output_dir = Path(output_dir).resolve()
    out_folder = output_dir / sys_id / f"{sys_id}_{lig_chain}"
    out_folder.mkdir(parents=True, exist_ok=True)

    cmd = [
        VINARDO,
        "--scoring", "2vinardo",
        "--receptor", receptor_file,
        "--ligand", ligand_file,
        "--out", str(out_folder),
        "--autobox",
        "--threads", str(threads),
    ]

    try:
        start = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600,
                                cwd=str(out_folder.parent))
        elapsed = time.time() - start

        log_txt = out_folder / "log.txt"
        with open(log_txt, "w") as fh:
            if result.stdout:
                fh.write("STDOUT:\n" + result.stdout)
            if result.stderr:
                fh.write("\nSTDERR:\n" + result.stderr)

        runtime_file = out_folder / "runtime.json"
        with open(runtime_file, "w") as fh:
            json.dump({"method": "vinardock", "threads": threads,
                       "returncode": result.returncode, "runtime_seconds": elapsed}, fh)

        if result.returncode != 0:
            return False

        log_file = out_folder / "log.csv"
        if log_file.exists() and log_file.stat().st_size > 50:
            return True
        output_files = list(out_folder.glob("*.pdbt")) + list(out_folder.glob("*.pdb"))
        if output_files and any(f.stat().st_size > 0 for f in output_files):
            return True
        return False
    except subprocess.TimeoutExpired:
        return False
    except Exception as e:
        print(f"Error: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(description="Run Vinardo on symmetry-corrected receptors")
    parser.add_argument("--receptor-dir", required=True)
    parser.add_argument("--ligand-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--system-list", required=True)
    parser.add_argument("--system-id", type=str, default=None)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    receptor_dir = Path(args.receptor_dir)
    ligand_dir = Path(args.ligand_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(args.system_list) as f:
        system_ids = [l.strip() for l in f if l.strip()]
    if args.system_id:
        system_ids = [s for s in system_ids if s == args.system_id]

    ok, fail, skip = 0, 0, 0
    for sys_id in system_ids:
        receptor_file = (receptor_dir / sys_id / f"{sys_id}_receptor.pdbt").resolve()
        if not receptor_file.exists():
            candidates = list((receptor_dir / sys_id).glob("*receptor*"))
            receptor_file = candidates[0] if candidates else None
        if not receptor_file:
            print(f"  No receptor: {sys_id}")
            fail += 1
            continue

        lig_dir = ligand_dir / sys_id
        if not lig_dir.exists():
            print(f"  No ligand dir: {sys_id}")
            fail += 1
            continue

        ligand_files = sorted(lig_dir.glob("*.pdbt"))
        if not ligand_files:
            print(f"  No ligand files: {sys_id}")
            fail += 1
            continue

        sys_ok = False
        for lf in ligand_files:
            lig_chain = lf.stem
            lf = lf.resolve()
            out_folder = output_dir / sys_id / f"{sys_id}_{lig_chain}"
            log_file = out_folder / "log.csv"

            if args.resume:
                log_csv = out_folder / "log.csv"
                if log_csv.exists() and log_csv.stat().st_size > 50:
                    skip += 1
                    sys_ok = True
                    continue
                output_files = list(out_folder.glob("*.pdbt")) + list(out_folder.glob("*.pdb"))
                if output_files and any(f.stat().st_size > 0 for f in output_files):
                    skip += 1
                    sys_ok = True
                    continue

            if run_vinardock(str(receptor_file), str(lf), str(output_dir),
                            sys_id, lig_chain, args.threads):
                ok += 1
                sys_ok = True
            else:
                print(f"  FAIL vinardo: {sys_id} {lig_chain}")
                fail += 1

        if not sys_ok:
            fail += 1

    print(f"Vinardo: OK={ok}, Skip={skip}, Fail={fail}")


if __name__ == "__main__":
    main()
