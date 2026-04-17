#!/usr/bin/env python
"""
Run the full AutoDock Vina pipeline on all single-ligand systems.

Steps:
1. Prepare receptor PDBQT files using `scripts/01_prepare_receptor_pdbqt_vina.py`
2. Prepare ligand PDBQT files using `scripts/02_prepare_ligand_pdbqt_vina.py`
3. Run `vina` using `scripts/03_run_vina.py`
4. Analyze results with `scripts/04_analyze_vina.py` and write predictions/posebusters CSVs
5. (optional) Plot Figure 1 replacing Boltz-1 with Vina using `scripts/05_plot_figure_with_vina.py`

Usage:
    python run_full_pipeline_vina.py [--step 1|2|3|4|5|all] [--threads N]
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm

SCRIPT_DIR = Path(__file__).parent
GT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth")
RECEPTOR_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/vina_inputs/receptors")
LIGAND_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/vina_inputs/ligands")
OUTPUT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/vina_outputs")
ANNOTATIONS = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/annotations.csv")
OBABEL = "obabel-25-07"
VINA = "vina"


def save_failed_list(failed: list, step_num: int):
    out_file = SCRIPT_DIR / f"failed_vina_step{step_num}.txt"
    with open(out_file, 'w') as f:
        for sid in failed:
            f.write(sid + '\n')
    if failed:
        print(f"\n  Failed systems ({len(failed)}) saved to: {out_file}")
    else:
        print(f"\n  No failures in step {step_num}!")


def step1_receptor(system_id: str) -> str:
    rec_dir = RECEPTOR_DIR / system_id
    rec_file = rec_dir / f"{system_id}_receptor.pdbqt"
    if rec_file.exists() and rec_file.stat().st_size > 1000:
        return "skip"
    cif = GT_DIR / system_id / "receptor.cif"
    if not cif.exists():
        return "no_cif"
    rec_dir.mkdir(parents=True, exist_ok=True)
    cmd = [OBABEL, str(cif), "-opdbqt", "-O", str(rec_file), "-h", "-xc", "-xr"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if r.returncode == 0 and rec_file.exists() and rec_file.stat().st_size > 1000:
            return "ok"
    except Exception:
        pass
    return "fail"


def step2_ligand(system_id: str) -> str:
    lig_dir = LIGAND_DIR / system_id
    if lig_dir.exists() and any(lig_dir.glob("*.pdbqt")):
        return "skip"
    sdf_dir = GT_DIR / system_id / "ligand_files"
    if not sdf_dir.exists():
        return "no_sdf"
    lig_dir.mkdir(parents=True, exist_ok=True)
    ok = 0
    for sdf_file in sorted(sdf_dir.glob("*.sdf")):
        chain = sdf_file.stem
        pdbqt = lig_dir / f"{chain}.pdbqt"
        cmd = [OBABEL, str(sdf_file), "-O", str(pdbqt), "-h"]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if r.returncode == 0 and pdbqt.exists() and pdbqt.stat().st_size > 0:
                ok += 1
        except Exception:
            pass
    return f"ok:{ok}" if ok > 0 else "fail"


def step3_run_vina(system_id: str, threads: int = 8) -> str:
    rec_files = list((RECEPTOR_DIR / system_id).glob("*receptor*.pdbqt"))
    if not rec_files:
        return "no_receptor"
    lig_files = list((LIGAND_DIR / system_id).glob("*.pdbqt"))
    if not lig_files:
        return "no_ligand"
    out_dir = OUTPUT_DIR / system_id
    ok = 0
    for lf in lig_files:
        chain = lf.stem
        out_folder = out_dir / f"{system_id}_{chain}"
        if out_folder.exists() and (out_folder / "log.txt").exists():
            if (out_folder / "log.txt").stat().st_size > 50:
                ok += 1
                continue
        out_folder.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, str(SCRIPT_DIR / "03_run_vina.py"),
               "--receptor-dir", str(RECEPTOR_DIR),
               "--ligand-dir", str(LIGAND_DIR),
               "--output-dir", str(OUTPUT_DIR),
               "--executable", VINA,
               "--system-id", system_id,
               "--threads", str(threads)]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
            # script creates per-ligand folders; check for any created folder
            if out_folder.exists() and any(out_folder.iterdir()):
                ok += 1
        except Exception:
            pass
    return f"ok:{ok}" if ok > 0 else "fail"


def step4_analyze(system_ids: list, resume: bool = False):
    cmd = [sys.executable, str(SCRIPT_DIR / "04_analyze_vina.py"), "--system-list", str(SCRIPT_DIR / "single_ligand_systems.txt")]
    if resume:
        cmd.append("--resume")
    r = subprocess.run(cmd)
    return r.returncode == 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--system-list", type=str,
                        default=str(SCRIPT_DIR / "single_ligand_systems.txt"))
    parser.add_argument("--start-index", type=int, default=0)
    parser.add_argument("--max-systems", type=int, default=None)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--step", type=str, choices=["1", "2", "3", "4", "5", "all"], default="all")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    with open(args.system_list) as f:
        system_ids = [l.strip() for l in f if l.strip()]
    system_ids = system_ids[args.start_index:]
    if args.max_systems:
        system_ids = system_ids[:args.max_systems]

    print(f"Processing {len(system_ids)} systems (threads={args.threads})")

    # Step 1
    if args.step in ("1", "all"):
        print("\n=== Step 1: Prepare receptors (Vina) ===")
        ok_list, fail_list = [], []
        for sid in tqdm(system_ids, desc="Receptors"):
            s = step1_receptor(sid)
            if s in ("ok", "skip"):
                ok_list.append(sid)
            else:
                fail_list.append(sid)
        print(f"  Done: {len(ok_list)} ok, {len(fail_list)} failed")
        save_failed_list(fail_list, 1)

    # Step 2
    if args.step in ("2", "all"):
        print("\n=== Step 2: Prepare ligands (Vina) ===")
        ok_list, fail_list = [], []
        for sid in tqdm(system_ids, desc="Ligands"):
            s = step2_ligand(sid)
            if s.startswith("ok") or s == "skip":
                ok_list.append(sid)
            else:
                fail_list.append(sid)
        print(f"  Done: {len(ok_list)} ok, {len(fail_list)} failed")
        save_failed_list(fail_list, 2)

    # Step 3
    if args.step in ("3", "all"):
        print("\n=== Step 3: Run Vina ===")
        ok_list, fail_list = [], []
        for sid in tqdm(system_ids, desc="Docking"):
            s = step3_run_vina(sid, args.threads)
            if s.startswith("ok") or s == "skip":
                ok_list.append(sid)
            else:
                fail_list.append(sid)
        print(f"  Done: {len(ok_list)} ok, {len(fail_list)} failed")
        save_failed_list(fail_list, 3)

    # Step 4
    if args.step in ("4", "all"):
        print("\n=== Step 4: Analyze results (Vina) ===")
        ok = step4_analyze(system_ids, resume=args.resume)
        if ok:
            # check predictions CSV
            output_csv = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/predictions/vina.csv")
            if output_csv.exists():
                results_df = pd.read_csv(output_csv)
                analyzed_systems = set(results_df['target'].unique()) if 'target' in results_df.columns else set()
                fail_list = [sid for sid in system_ids if sid not in analyzed_systems]
                ok_list = [sid for sid in system_ids if sid in analyzed_systems]
                print(f"\n  Done: {len(ok_list)} ok, {len(fail_list)} failed")
                save_failed_list(fail_list, 4)
            else:
                print("  Analysis failed - no output CSV created")
                save_failed_list(system_ids, 4)
        else:
            print("  Analysis script exited with non-zero code")
            save_failed_list(system_ids, 4)

    # Step 5: plotting
    if args.step in ("5", "all"):
        print("\n=== Step 5: Plot figures (Vina) ===")
        env = os.environ.copy()
        env["MPLBACKEND"] = "Agg"
        cmd = [
            sys.executable, "-m", "conda", "run", "-n", "runs_n_poses",
            "python", str(SCRIPT_DIR / "05_plot_figure_with_vina.py")
        ]
        r = subprocess.run(cmd, env=env)
        if r.returncode == 0:
            print("\n  Done: figures generated successfully")
            save_failed_list([], 5)
        else:
            print(f"\n  Failed: exit code {r.returncode}")
            save_failed_list(system_ids, 5)


if __name__ == "__main__":
    main()
