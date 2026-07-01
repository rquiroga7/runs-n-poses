#!/usr/bin/env python
"""
Run the full vinardo pipeline on all single-ligand systems.

Steps:
1. Prepare receptor PDBT files (obabel -h -xc -xr)
2. Prepare ligand PDBT files 
3. Run vinardock-26-04 --scoring 2vinardo --autobox --threads N
4. Analyze results with ost and output predictions CSV
5. Plot Figure 1E with vinardo results added

Each step outputs a list of failed accessions to scripts/failed_step{N}.txt

Usage:
    python run_full_pipeline.py [--step 1|2|3|4|5|all] [--threads N]
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
GT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth")
RECEPTOR_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_inputs/receptors")
LIGAND_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_inputs/ligands")
OUTPUT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_outputs")
ANNOTATIONS = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/annotations.csv")
OBABEL = "obabel-25-07"
VINARDOCK = "vinardock-26-04"
RUNS_N_POSES_PYTHON = "/home/rquiroga/anaconda3/envs/runs_n_poses/bin/python"


def _bootstrap_project_python() -> None:
    """Re-exec into the project interpreter if the current python lacks deps."""
    if os.environ.get("RUNS_N_POSES_BOOTSTRAPPED") == "1":
        return
    try:
        import pandas  # noqa: F401
        return
    except ModuleNotFoundError:
        pass

    if sys.executable == RUNS_N_POSES_PYTHON:
        raise

    os.environ["RUNS_N_POSES_BOOTSTRAPPED"] = "1"
    os.execv(RUNS_N_POSES_PYTHON, [RUNS_N_POSES_PYTHON, str(Path(__file__).resolve()), *sys.argv[1:]])


_bootstrap_project_python()

import pandas as pd
from tqdm import tqdm
from pandas.errors import EmptyDataError


def save_failed_list(failed: list, step_num: int):
    """Save list of failed accessions to file."""
    out_file = SCRIPT_DIR / f"failed_step{step_num}.txt"
    with open(out_file, 'w') as f:
        for sid in failed:
            f.write(sid + '\n')
    if failed:
        print(f"\n  Failed systems ({len(failed)}) saved to: {out_file}")
    else:
        print(f"\n  No failures in step {step_num}!")


def step1_receptor(system_id: str) -> str:
    """Prepare receptor PDBT."""
    rec_dir = RECEPTOR_DIR / system_id
    rec_file = rec_dir / f"{system_id}_receptor.pdbt"
    if rec_file.exists() and rec_file.stat().st_size > 1000:
        return "skip"
    cif = GT_DIR / system_id / "receptor.cif"
    if not cif.exists():
        return "no_cif"
    rec_dir.mkdir(parents=True, exist_ok=True)
    r = subprocess.run([OBABEL, str(cif), "-opdbt", "-O", str(rec_file), "-h", "-xc", "-xr"],
                       capture_output=True, text=True, timeout=120)
    if r.returncode == 0 and rec_file.exists() and rec_file.stat().st_size > 1000:
        return "ok"
    return "fail"


def step2_ligand(system_id: str) -> str:
    """Prepare ligand PDBT files from ground truth SDFs."""
    lig_dir = LIGAND_DIR / system_id
    if lig_dir.exists() and any(lig_dir.glob("*.pdbt")):
        return "skip"
    sdf_dir = GT_DIR / system_id / "ligand_files"
    if not sdf_dir.exists():
        return "no_sdf"
    lig_dir.mkdir(parents=True, exist_ok=True)
    ok = 0
    for sdf_file in sorted(sdf_dir.glob("*.sdf")):
        chain = sdf_file.stem
        pdbt = lig_dir / f"{chain}.pdbt"
        cmd = [OBABEL, str(sdf_file), "-O", str(pdbt), "-h"]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if r.returncode == 0 and pdbt.exists() and pdbt.stat().st_size > 0:
                ok += 1
        except Exception:
            pass
    return f"ok:{ok}" if ok > 0 else "fail"


def step3_dock(system_id: str, threads: int = 10) -> str:
    """Run vinardock docking."""
    rec_files = list((RECEPTOR_DIR / system_id).glob("*receptor*.pdbt"))
    if not rec_files:
        return "no_receptor"
    lig_files = list((LIGAND_DIR / system_id).glob("*.pdbt"))
    if not lig_files:
        return "no_ligand"
    out_dir = OUTPUT_DIR / system_id
    ok = 0
    for lf in lig_files:
        chain = lf.stem
        out_folder = out_dir / f"{system_id}_{chain}"
        if out_folder.exists() and (out_folder / "log.csv").exists():
            if (out_folder / "log.csv").stat().st_size > 50:
                ok += 1
                continue
        out_folder.mkdir(parents=True, exist_ok=True)
        cmd = [VINARDOCK, "--scoring", "2vinardo", "--receptor", str(rec_files[0]),
               "--ligand", str(lf), "--out", str(out_folder), "--autobox", "--threads", str(threads)]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            log = out_folder / "log.csv"
            if r.returncode == 0 and log.exists() and log.stat().st_size > 50:
                ok += 1
        except Exception:
            pass
    return f"ok:{ok}" if ok > 0 else "fail"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--system-list", type=str,
                        default="/home/rquiroga/github/runs-n-poses/scripts/single_ligand_systems.txt")
    parser.add_argument("--start-index", type=int, default=0)
    parser.add_argument("--max-systems", type=int, default=None)
    parser.add_argument("--threads", type=int, default=10)
    parser.add_argument("--step", type=str, choices=["1", "2", "3", "4", "5", "all"], default="all")
    parser.add_argument("--resume", action="store_true",
                        help="When running step 4, resume from existing analysis and skip already processed systems")
    args = parser.parse_args()

    with open(args.system_list) as f:
        system_ids = [l.strip() for l in f if l.strip()]
    system_ids = system_ids[args.start_index:]
    if args.max_systems:
        system_ids = system_ids[:args.max_systems]

    print(f"Processing {len(system_ids)} systems (threads={args.threads})")

    # Step 1: Prepare receptors (delegate to external script)
    if args.step in ("1", "all"):
        print("\n=== Step 1: Prepare receptors ===")
        ok_list, fail_list = [], []
        for sid in tqdm(system_ids, desc="Receptors"):
            cmd = [RUNS_N_POSES_PYTHON, str(SCRIPT_DIR / "01_prepare_receptor_pdbqt.py"),
                   "--ground-truth-dir", str(GT_DIR),
                   "--output-dir", str(RECEPTOR_DIR),
                   "--system-id", sid]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode == 0:
                ok_list.append(sid)
            else:
                fail_list.append(sid)
        print(f"  Done: {len(ok_list)} ok, {len(fail_list)} failed")
        save_failed_list(fail_list, 1)

    # Step 2: Prepare ligands (delegate to external script)
    if args.step in ("2", "all"):
        print("\n=== Step 2: Prepare ligands ===")
        ok_list, fail_list = [], []
        for sid in tqdm(system_ids, desc="Ligands"):
            cmd = [RUNS_N_POSES_PYTHON, str(SCRIPT_DIR / "02_prepare_ligand_pdbqt.py"),
                   "--ground-truth-dir", str(GT_DIR),
                   "--output-dir", str(LIGAND_DIR),
                   "--system-id", sid]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode == 0:
                ok_list.append(sid)
            else:
                fail_list.append(sid)
        print(f"  Done: {len(ok_list)} ok, {len(fail_list)} failed")
        save_failed_list(fail_list, 2)

    # Step 3: Run vinardock (delegate to external script)
    if args.step in ("3", "all"):
        print("\n=== Step 3: Run vinardock ===")
        ok_list, fail_list = [], []
        for sid in tqdm(system_ids, desc="Docking"):
            # If resume requested, skip systems where all ligand runs already have a
            # non-empty log.csv (heuristic for completed docking).
            if args.resume:
                lig_dir = LIGAND_DIR / sid
                if lig_dir.exists():
                    lig_files = list(lig_dir.glob("*.pdbt"))
                else:
                    lig_files = []
                if lig_files:
                    all_done = True
                    for lf in lig_files:
                        out_folder = OUTPUT_DIR / sid / f"{sid}_{lf.stem}"
                        log_file = out_folder / "log.csv"
                        if not (out_folder.exists() and log_file.exists() and log_file.stat().st_size > 50):
                            all_done = False
                            break
                    if all_done:
                        tqdm.write(f"Skipping {sid} (already docked)")
                        ok_list.append(sid)
                        continue
                        cmd = [RUNS_N_POSES_PYTHON, str(SCRIPT_DIR / "03_run_vinardo.py"),
                   "--receptor-dir", str(RECEPTOR_DIR),
                   "--ligand-dir", str(LIGAND_DIR),
                   "--output-dir", str(OUTPUT_DIR),
                   "--system-id", sid,
                   "--threads", str(args.threads)]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode == 0:
                ok_list.append(sid)
            else:
                fail_list.append(sid)
        print(f"  Done: {len(ok_list)} ok, {len(fail_list)} failed")
        save_failed_list(fail_list, 3)

    # Step 4: Analyze results
    if args.step in ("4", "all"):
        print("\n=== Step 4: Analyze results ===")
        # Run the analysis script (let output pass through for progress bars)
        cmd = [RUNS_N_POSES_PYTHON, str(SCRIPT_DIR / "04_analyze_vinardo.py"),
               "--system-list", args.system_list]
        if args.start_index > 0:
            cmd += ["--start-index", str(args.start_index)]
        if args.max_systems:
            cmd += ["--max-systems", str(args.max_systems)]
        if args.resume:
            cmd += ["--resume"]

        # Run with output passing through to terminal (so progress bars show)
        r = subprocess.run(cmd)
        
        # Read the results CSV to determine success/failure
        output_csv = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/predictions/vinardock_2vinardo.csv")
        if output_csv.exists():
            try:
                results_df = pd.read_csv(output_csv)
            except EmptyDataError:
                results_df = None
            if results_df is not None and 'target' in results_df.columns:
                analyzed_systems = set(results_df['target'].unique())
                fail_list = [sid for sid in system_ids if sid not in analyzed_systems]
                ok_list = [sid for sid in system_ids if sid in analyzed_systems]
                print(f"\n  Done: {len(ok_list)} ok, {len(fail_list)} failed")
                save_failed_list(fail_list, 4)
            else:
                print("  Analysis produced no parseable output CSV")
                save_failed_list(system_ids, 4)
        else:
            print("  Analysis failed - no output CSV created")
            save_failed_list(system_ids, 4)

    # Step 5: Plot figures
    if args.step in ("5", "all"):
        print("\n=== Step 5: Plot figures ===")
        # plotting.py needs cmap/ost from the current environment, and Agg backend for headless
        env = os.environ.copy()
        env["MPLBACKEND"] = "Agg"
        cmd = [RUNS_N_POSES_PYTHON, str(SCRIPT_DIR / "05_plot_figure_with_vinardo.py")]
        r = subprocess.run(cmd, env=env)
        if r.returncode == 0:
            print(f"\n  Done: figures generated successfully")
            save_failed_list([], 5)
        else:
            print(f"\n  Failed: exit code {r.returncode}")
            save_failed_list(system_ids, 5)


if __name__ == "__main__":
    main()
