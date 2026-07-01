#!/usr/bin/env python
"""
Step 3: Run vinardock-26-04 --scoring 2vinardo for each complex.

This script runs vinardock-26-04 --scoring 2vinardo --autobox for each system,
using the prepared receptor and ligand PDBT files.

Usage:
    python 03_run_vinardo.py --receptor-dir /path/to/receptors \
                             --ligand-dir /path/to/ligands \
                             --output-dir /path/to/outputs
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path
import time
import json

import pandas as pd
from tqdm import tqdm


def run_vinardock(receptor_file: str, ligand_file: str, output_dir: str,
                  system_id: str, ligand_chain: str,
                  executable: str = "vinardock-26-04",
                  threads: int = 4) -> bool:
    """
    Run vinardock-26-04 --scoring 2vinardo --autobox for a single receptor-ligand pair.
    """
    system_output_dir = Path(output_dir) / system_id
    system_output_dir.mkdir(parents=True, exist_ok=True)

    out_folder = system_output_dir / f"{system_id}_{ligand_chain}"
    out_folder.mkdir(exist_ok=True)

    cmd = [
        executable,
        "--scoring", "2vinardo",
        "--receptor", receptor_file,
        "--ligand", ligand_file,
        "--out", str(out_folder),
        "--autobox",
        "--threads", str(threads)
    ]

    try:
        start_time = time.time()
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,
            cwd=str(system_output_dir)
        )
        end_time = time.time()
        elapsed = end_time - start_time

        # Save stdout/stderr to a log file for diagnostics
        log_txt = out_folder / "log.txt"
        try:
            with open(log_txt, 'w') as fh:
                if result.stdout:
                    fh.write('STDOUT:\n')
                    fh.write(result.stdout)
                if result.stderr:
                    fh.write('\nSTDERR:\n')
                    fh.write(result.stderr)
        except Exception:
            pass

        # Record runtime information for later analysis
        try:
            run_meta = {
                "method": "vinardock",
                "threads": int(threads),
                "returncode": result.returncode,
                "runtime_seconds": float(elapsed),
            }
            runtime_file = out_folder / "runtime.json"
            with open(runtime_file, 'w') as fh:
                json.dump(run_meta, fh)
        except Exception:
            pass

        if result.returncode != 0:
            print(f"Warning: vinardock failed for {system_id} {ligand_chain}:")
            print(f"  stderr: {result.stderr[:500]}")
            return False

        # Check if output files were created
        output_files = list(out_folder.glob("*.pdbt")) + list(out_folder.glob("*.pdb"))
        # Also check if log.csv has content
        log_file = out_folder / "log.csv"
        if log_file.exists() and log_file.stat().st_size > 50:
            return True
        if output_files and any(f.stat().st_size > 0 for f in output_files):
            return True

        print(f"Warning: No output files created for {system_id} {ligand_chain}")
        return False

    except subprocess.TimeoutExpired:
        print(f"Warning: vinardock timed out for {system_id} {ligand_chain}")
        return False
    except Exception as e:
        print(f"Error running vinardock for {system_id} {ligand_chain}: {e}")
        return False


def process_system(system_id: str, receptor_dir: Path, ligand_dir: Path,
                   output_dir: Path, executable: str,
                   threads: int = 4) -> bool:
    """Process a single system by running vinardock for all its ligands."""
    receptor_file = receptor_dir / system_id / f"{system_id}_receptor.pdbt"
    if not receptor_file.exists():
        candidates = list((receptor_dir / system_id).glob("*receptor*"))
        if candidates:
            receptor_file = candidates[0]
        else:
            print(f"Warning: Receptor file not found for {system_id}")
            return False

    ligand_dir_system = ligand_dir / system_id
    if not ligand_dir_system.exists():
        print(f"Warning: Ligand directory not found for {system_id}")
        return False

    ligand_files = list(ligand_dir_system.glob("*.pdbt"))
    if not ligand_files:
        print(f"Warning: No ligand files found for {system_id}")
        return False

    success = False
    for ligand_file in ligand_files:
        ligand_chain = ligand_file.stem
        print(f"Running vinardock for {system_id}, ligand {ligand_chain}...")

        if run_vinardock(
            str(receptor_file),
            str(ligand_file),
            str(output_dir),
            system_id,
            ligand_chain,
            executable,
            threads
        ):
            success = True
        else:
            print(f"  Failed for ligand {ligand_chain}")

    return success


def main():
    parser = argparse.ArgumentParser(description="Run vinardock-26-04 --scoring 2vinardo for each complex")
    parser.add_argument(
        "--receptor-dir",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_inputs/receptors",
        help="Path to directory containing receptor PDBT files"
    )
    parser.add_argument(
        "--ligand-dir",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_inputs/ligands",
        help="Path to directory containing ligand PDBT files"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_outputs",
        help="Path to output directory for docking results"
    )
    parser.add_argument(
        "--executable",
        type=str,
        default="vinardock-26-04",
        help="Path to vinardock executable"
    )
    parser.add_argument(
        "--system-id",
        type=str,
        default=None,
        help="Process only a specific system ID (optional)"
    )
    parser.add_argument(
        "--annotations",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/annotations.csv",
        help="Path to annotations.csv file"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads per docking job"
    )

    args = parser.parse_args()

    receptor_dir = Path(args.receptor_dir)
    ligand_dir = Path(args.ligand_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if os.path.exists(args.annotations):
        annotations = pd.read_csv(args.annotations)
        if args.system_id:
            system_ids = [args.system_id]
        else:
            system_ids = sorted(annotations["system_id"].unique())
    else:
        system_ids = sorted([d.name for d in receptor_dir.iterdir() if d.is_dir()])
        if args.system_id:
            system_ids = [s for s in system_ids if s == args.system_id]

    print(f"Processing {len(system_ids)} systems...")

    success_count = 0
    fail_count = 0

    for system_id in tqdm(system_ids, desc="Running vinardock"):
        if process_system(
            system_id,
            receptor_dir,
            ligand_dir,
            output_dir,
            args.executable,
            args.threads
        ):
            success_count += 1
        else:
            fail_count += 1

    print(f"\nDone! Success: {success_count}, Failed: {fail_count}")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()
