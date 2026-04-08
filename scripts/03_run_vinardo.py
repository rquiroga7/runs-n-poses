#!/usr/bin/env python
"""
Step 3: Run 2vinardo-mar5_autobox for each complex.

This script runs the 2vinardo-mar5_autobox executable for each system,
using the prepared receptor and ligand PDBQT files.

Usage:
    python 03_run_vinardo.py --receptor-dir /path/to/receptors \
                             --ligand-dir /path/to/ligands \
                             --output-dir /path/to/outputs \
                             --config /path/to/config.fijo
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm


def run_vinardo(receptor_file: str, ligand_file: str, output_dir: str, 
                config_file: str, system_id: str, ligand_chain: str,
                executable: str = "2vinardo-mar5_autobox") -> bool:
    """
    Run 2vinardo-mar5_autobox for a single receptor-ligand pair.
    
    Args:
        receptor_file: Path to receptor PDBQT file
        ligand_file: Path to ligand PDBQT file
        output_dir: Directory to save output files
        config_file: Path to config.fijo file
        system_id: System ID for naming output
        ligand_chain: Ligand chain identifier
        executable: Path to 2vinardo-mar5_autobox executable
    """
    # Create output subdirectory for this system
    system_output_dir = Path(output_dir) / system_id
    system_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Output file name
    output_file = system_output_dir / f"{system_id}_{ligand_chain}_output.pdbqt"
    
    cmd = [
        executable,
        "--receptor", receptor_file,
        "--ligand", ligand_file,
        "--config", config_file,
        "--out", str(output_file)
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,  # 10 minutes timeout
            cwd=str(system_output_dir)
        )
        
        if result.returncode != 0:
            print(f"Warning: 2vinardo-mar5 failed for {system_id} {ligand_chain}:")
            print(f"  stderr: {result.stderr[:500]}")
            return False
        
        # Check if output file was created
        if not output_file.exists() or output_file.stat().st_size == 0:
            print(f"Warning: No output file created for {system_id} {ligand_chain}")
            return False
        
        return True
        
    except subprocess.TimeoutExpired:
        print(f"Warning: 2vinardo-mar5 timed out for {system_id} {ligand_chain}")
        return False
    except Exception as e:
        print(f"Error running 2vinardo-mar5 for {system_id} {ligand_chain}: {e}")
        return False


def process_system(system_id: str, receptor_dir: Path, ligand_dir: Path,
                   output_dir: Path, config_file: str, executable: str) -> bool:
    """
    Process a single system by running 2vinardo-mar5 for all its ligands.
    """
    receptor_file = receptor_dir / system_id / f"{system_id}_receptor.pdbqt"
    ligand_dir_system = ligand_dir / system_id
    
    if not receptor_file.exists():
        print(f"Warning: Receptor file not found for {system_id}")
        return False
    
    if not ligand_dir_system.exists():
        print(f"Warning: Ligand directory not found for {system_id}")
        return False
    
    # Get all ligand PDBQT files for this system
    ligand_files = list(ligand_dir_system.glob("*.pdbqt"))
    
    if not ligand_files:
        print(f"Warning: No ligand files found for {system_id}")
        return False
    
    success = False
    
    for ligand_file in ligand_files:
        ligand_chain = ligand_file.stem  # e.g., "1.B"
        
        print(f"Running 2vinardo-mar5 for {system_id}, ligand {ligand_chain}...")
        
        if run_vinardo(
            str(receptor_file),
            str(ligand_file),
            str(output_dir),
            config_file,
            system_id,
            ligand_chain,
            executable
        ):
            success = True
        else:
            print(f"  Failed for ligand {ligand_chain}")
    
    return success


def main():
    parser = argparse.ArgumentParser(description="Run 2vinardo-mar5_autobox for each complex")
    parser.add_argument(
        "--receptor-dir",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_inputs/receptors",
        help="Path to directory containing receptor PDBQT files"
    )
    parser.add_argument(
        "--ligand-dir",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_inputs/ligands",
        help="Path to directory containing ligand PDBQT files"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_outputs",
        help="Path to output directory for docking results"
    )
    parser.add_argument(
        "--config",
        type=str,
        default="/home/rquiroga/github/runs-n-poses/config.fijo",
        help="Path to config.fijo file"
    )
    parser.add_argument(
        "--executable",
        type=str,
        default="2vinardo-mar5_autobox",
        help="Path to 2vinardo-mar5_autobox executable"
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
        default="/home/rquiroga/Datasets/runs-n-Poses/annotations.csv",
        help="Path to annotations.csv file (to determine which ligands to process)"
    )
    
    args = parser.parse_args()
    
    receptor_dir = Path(args.receptor_dir)
    ligand_dir = Path(args.ligand_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load annotations to get list of systems
    if os.path.exists(args.annotations):
        annotations = pd.read_csv(args.annotations)
        if args.system_id:
            system_ids = [args.system_id]
        else:
            system_ids = sorted(annotations["system_id"].unique())
    else:
        # Use directories instead
        system_ids = sorted([d.name for d in receptor_dir.iterdir() if d.is_dir()])
        if args.system_id:
            system_ids = [s for s in system_ids if s == args.system_id]
    
    print(f"Processing {len(system_ids)} systems...")
    
    success_count = 0
    fail_count = 0
    
    for system_id in tqdm(system_ids, desc="Running 2vinardo-mar5"):
        if process_system(
            system_id,
            receptor_dir,
            ligand_dir,
            output_dir,
            args.config,
            args.executable
        ):
            success_count += 1
        else:
            fail_count += 1
    
    print(f"\nDone! Success: {success_count}, Failed: {fail_count}")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()
