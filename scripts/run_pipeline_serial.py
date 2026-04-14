#!/usr/bin/env python
"""
Run steps 1-3 of the vinardo pipeline serially (one system at a time).

Vinardock handles its own parallelization via --threads flag.

Usage:
    python run_pipeline_serial.py [--start-index N] [--max-systems N] [--threads N]
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm

GT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth")
RECEPTOR_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_inputs/receptors")
LIGAND_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_inputs/ligands")
OUTPUT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_outputs")
ANNOTATIONS = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/annotations.csv")
OBABEL = "obabel-25-07"
VINARDOCK = "/home/rquiroga/.local/bin/vinardock-26-mar"


def run_step1(system_id: str) -> str:
    """Prepare receptor PDBT. Returns status."""
    receptor_dir = RECEPTOR_DIR / system_id
    receptor_file = receptor_dir / f"{system_id}_receptor.pdbt"
    if receptor_file.exists() and receptor_file.stat().st_size > 1000:
        return "skip"

    cif = GT_DIR / system_id / "receptor.cif"
    if not cif.exists():
        return "no_cif"

    receptor_dir.mkdir(parents=True, exist_ok=True)
    cmd = [OBABEL, str(cif), "-opdbt", "-O", str(receptor_file), "-d", "-xc", "-xr"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if r.returncode == 0 and receptor_file.exists() and receptor_file.stat().st_size > 1000:
            return "ok"
        return "fail"
    except Exception as e:
        return f"error:{e}"


def run_step2(system_id: str) -> str:
    """Prepare ligand PDBT files. Returns status."""
    ligand_dir = LIGAND_DIR / system_id
    if ligand_dir.exists() and any(ligand_dir.glob("*.pdbt")):
        return "skip"

    annot_df = pd.read_csv(ANNOTATIONS)
    rows = annot_df[annot_df["system_id"] == system_id]
    if rows.empty:
        return "no_annot"

    ligand_dir.mkdir(parents=True, exist_ok=True)
    ok_count = 0
    for _, row in rows.iterrows():
        chain = row["ligand_instance_chain"]
        smiles = row.get("ligand_smiles", "")
        pdbt_file = ligand_dir / f"{chain}.pdbt"

        if pd.notna(smiles) and smiles:
            cmd = [OBABEL, "-ismi", "-:" + smiles, "-opdbt", "-O", str(pdbt_file), "--gen3d", "-d"]
        else:
            sdf = GT_DIR / system_id / "ligand_files" / f"{chain}.sdf"
            if not sdf.exists():
                continue
            cmd = [OBABEL, str(sdf), "-O", str(pdbt_file)]

        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if r.returncode == 0 and pdbt_file.exists() and pdbt_file.stat().st_size > 0:
                ok_count += 1
        except Exception:
            pass

    return f"ok:{ok_count}" if ok_count > 0 else "fail"


def run_step3(system_id: str, threads: int = 10) -> str:
    """Run vinardock docking. Returns status."""
    receptor_files = list((RECEPTOR_DIR / system_id).glob("*receptor*.pdbt"))
    if not receptor_files:
        return "no_receptor"
    receptor_file = receptor_files[0]

    ligand_files = list((LIGAND_DIR / system_id).glob("*.pdbt"))
    if not ligand_files:
        return "no_ligand"

    out_system_dir = OUTPUT_DIR / system_id
    ok_count = 0
    for ligand_file in ligand_files:
        chain = ligand_file.stem
        out_folder = out_system_dir / f"{system_id}_{chain}"
        if out_folder.exists() and (out_folder / "log.csv").exists():
            if (out_folder / "log.csv").stat().st_size > 50:
                ok_count += 1
                continue

        out_folder.mkdir(parents=True, exist_ok=True)
        cmd = [VINARDOCK, "--scoring", "2vinardo", "--receptor", str(receptor_file),
               "--ligand", str(ligand_file), "--out", str(out_folder),
               "--autobox", "--threads", str(threads)]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            log = out_folder / "log.csv"
            if r.returncode == 0 and log.exists() and log.stat().st_size > 50:
                ok_count += 1
        except Exception:
            pass

    return f"ok:{ok_count}" if ok_count > 0 else "fail"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--system-list", type=str,
                        default="/home/rquiroga/github/runs-n-poses/scripts/single_ligand_systems.txt")
    parser.add_argument("--start-index", type=int, default=0, help="Start from this index")
    parser.add_argument("--max-systems", type=int, default=None, help="Max systems to process")
    parser.add_argument("--threads", type=int, default=10, help="Threads per vinardock job")
    parser.add_argument("--step", type=str, choices=["1", "2", "3", "all"], default="all")
    args = parser.parse_args()

    with open(args.system_list) as f:
        system_ids = [line.strip() for line in f if line.strip()]

    system_ids = system_ids[args.start_index:]
    if args.max_systems:
        system_ids = system_ids[:args.max_systems]

    print(f"Processing {len(system_ids)} systems (threads={args.threads})")

    if args.step in ("1", "all"):
        print("\n=== Step 1: Prepare receptors ===")
        ok1, fail1 = 0, 0
        for sid in tqdm(system_ids, desc="Receptors"):
            status = run_step1(sid)
            if status.startswith("ok") or status == "skip":
                ok1 += 1
            else:
                fail1 += 1
        print(f"Step 1 done: {ok1} ok, {fail1} failed")

    if args.step in ("2", "all"):
        print("\n=== Step 2: Prepare ligands ===")
        ok2, fail2 = 0, 0
        for sid in tqdm(system_ids, desc="Ligands"):
            status = run_step2(sid)
            if status.startswith("ok") or status == "skip":
                ok2 += 1
            else:
                fail2 += 1
        print(f"Step 2 done: {ok2} ok, {fail2} failed")

    if args.step in ("3", "all"):
        print("\n=== Step 3: Run vinardock ===")
        ok3, fail3 = 0, 0
        for sid in tqdm(system_ids, desc="Docking"):
            status = run_step3(sid, args.threads)
            if status.startswith("ok") or status == "skip":
                ok3 += 1
            else:
                fail3 += 1
        print(f"Step 3 done: {ok3} ok, {fail3} failed")


if __name__ == "__main__":
    main()
