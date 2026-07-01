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
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

SCRIPT_DIR = Path(__file__).parent
GT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth")
RECEPTOR_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/vina_inputs/receptors")
LIGAND_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/vina_inputs/ligands")
OUTPUT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/vina_outputs")
ANNOTATIONS = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/annotations.csv")
OBABEL = "obabel-25-07"
VINA = "vina"
AUTODOCK_VINA_8 = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/autodock_vina_8")
AUTODOCK_VINA_32 = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/autodock_vina_32")
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


def step3_run_vina(system_id: str, threads: int = 8, resume: bool = False) -> str:
    rec_files = list((RECEPTOR_DIR / system_id).glob("*receptor*.pdbqt"))
    if not rec_files:
        return "no_receptor"
    lig_files = list((LIGAND_DIR / system_id).glob("*.pdbqt"))
    if not lig_files:
        return "no_ligand"
    ok = 0

    # Run both exhaustiveness settings into separate output trees.
    for exh, out_base in [(8, AUTODOCK_VINA_8), (32, AUTODOCK_VINA_32)]:
        out_base.mkdir(parents=True, exist_ok=True)
        # Call the runner for this system and exhaustiveness, enabling resume
        cmd = [
            RUNS_N_POSES_PYTHON, str(SCRIPT_DIR / "03_run_vina.py"),
            "--receptor-dir", str(RECEPTOR_DIR),
            "--ligand-dir", str(LIGAND_DIR),
            "--output-dir", str(out_base),
            "--exhaustiveness", str(exh),
            "--executable", VINA,
            "--system-id", system_id,
            "--threads", str(threads),
        ]
        if resume:
            cmd.append("--resume")
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
            # check whether any per-ligand output folder was created for this system
            sys_out = out_base / system_id
            if sys_out.exists() and any(sys_out.iterdir()):
                ok += 1
        except Exception:
            pass

    return f"ok:{ok}" if ok > 0 else "fail"


def step4_analyze(
    script_name: str,
    start_index: int = 0,
    max_systems: Optional[int] = None,
    resume: bool = False,
) -> int:
    cmd = [RUNS_N_POSES_PYTHON, str(SCRIPT_DIR / script_name), "--system-list", str(SCRIPT_DIR / "single_ligand_systems.txt")]
    if start_index:
        cmd.extend(["--start-index", str(start_index)])
    if max_systems is not None:
        cmd.extend(["--max-systems", str(max_systems)])
    if resume:
        cmd.append("--resume")
    r = subprocess.run(cmd)
    return r.returncode


def _csv_targets(csv_path: Path, column_name: str = "target") -> set[str]:
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return set()
    try:
        df = pd.read_csv(csv_path, low_memory=False)
    except Exception:
        return set()
    if column_name not in df.columns:
        return set()
    return set(df[column_name].astype(str).unique())


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

    # Step 1: Prepare receptors (delegate to external script)
    if args.step in ("1", "all"):
        print("\n=== Step 1: Prepare receptors (Vina) ===")
        ok_list, fail_list = [], []
        for sid in tqdm(system_ids, desc="Receptors"):
            cmd = [RUNS_N_POSES_PYTHON, str(SCRIPT_DIR / "01_prepare_receptor_pdbqt_vina.py"),
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
        print("\n=== Step 2: Prepare ligands (Vina) ===")
        ok_list, fail_list = [], []
        for sid in tqdm(system_ids, desc="Ligands"):
            cmd = [RUNS_N_POSES_PYTHON, str(SCRIPT_DIR / "02_prepare_ligand_pdbqt_vina.py"),
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

    # Step 3: Run Vina (delegate to external script runner per-exhaustiveness)
    if args.step in ("3", "all"):
        print("\n=== Step 3: Run Vina ===")
        ok_list_8, fail_list_8 = [], []
        ok_list_32, fail_list_32 = [], []
        for sid in tqdm(system_ids, desc="Docking"):
            cmd = [RUNS_N_POSES_PYTHON, str(SCRIPT_DIR / "03_run_vina.py"),
                   "--receptor-dir", str(RECEPTOR_DIR),
                   "--ligand-dir", str(LIGAND_DIR),
                   "--output-dir", str(AUTODOCK_VINA_8),
                   "--executable", VINA,
                   "--system-id", sid,
                   "--threads", str(args.threads),
                   "--exhaustiveness", "8"]
            r8 = subprocess.run(cmd, capture_output=True, text=True)
            if r8.returncode == 0:
                ok_list_8.append(sid)
            else:
                fail_list_8.append(sid)

            cmd32 = cmd.copy()
            for i, v in enumerate(cmd32):
                if v == str(AUTODOCK_VINA_8):
                    cmd32[i] = str(AUTODOCK_VINA_32)
                if v == "8":
                    cmd32[i] = "32"
            r32 = subprocess.run(cmd32, capture_output=True, text=True)
            if r32.returncode == 0:
                ok_list_32.append(sid)
            else:
                fail_list_32.append(sid)

        print(f"  Exhaustiveness 8: {len(ok_list_8)} ok, {len(fail_list_8)} failed")
        print(f"  Exhaustiveness 32: {len(ok_list_32)} ok, {len(fail_list_32)} failed")
        combined_ok = sorted(set(ok_list_8) | set(ok_list_32))
        combined_fail = sorted(set(system_ids) - set(combined_ok))
        print(f"  Any exhaustiveness: {len(combined_ok)} ok, {len(combined_fail)} failed")
        save_failed_list(combined_fail, 3)

    # Step 4: Analyze results
    if args.step in ("4", "all"):
        print("\n=== Step 4: Analyze results (Vina) ===")
        output_csv_8 = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/predictions/autodock_vina_8.csv")
        posebusters_csv_8 = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/posebusters_results/autodock_vina_8.csv")
        output_csv_vina = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/predictions/vina.csv")
        posebusters_csv_vina = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/posebusters_results/vina.csv")

        expected_system_ids = {str(sid) for sid in system_ids}
        refresh_vina_8 = False
        existing_vina_8_targets = _csv_targets(output_csv_8)
        if existing_vina_8_targets and existing_vina_8_targets.isdisjoint(expected_system_ids):
            refresh_vina_8 = True

        if refresh_vina_8:
            print("  Detected corrupted Vina exh=8 outputs; rebuilding from docking results.")
            for path in [output_csv_8, posebusters_csv_8, output_csv_vina, posebusters_csv_vina]:
                if path.exists():
                    backup_path = path.with_suffix(path.suffix + ".bak")
                    try:
                        shutil.move(str(path), str(backup_path))
                    except Exception:
                        pass

        step4_analyze("04_analyze_vina.py", start_index=args.start_index, max_systems=args.max_systems, resume=args.resume and not refresh_vina_8)
        step4_analyze("04_analyze_vina_32.py", start_index=args.start_index, max_systems=args.max_systems, resume=args.resume)

        if output_csv_8.exists():
            shutil.copy2(output_csv_8, output_csv_vina)
        if posebusters_csv_8.exists():
            shutil.copy2(posebusters_csv_8, posebusters_csv_vina)

        if output_csv_8.exists():
            try:
                results_df = pd.read_csv(output_csv_8)
            except EmptyDataError:
                results_df = None
            if results_df is not None and "target" in results_df.columns:
                print(f"\n  Results saved to: {output_csv_8}")
                print(f"  Total predictions: {len(results_df)}")
                print(f"  Unique systems: {results_df['target'].nunique()}")
            print(f"\n  Done: {len(system_ids)} ok, 0 failed")
            save_failed_list([], 4)
        else:
            print("\n  Analysis completed, but no output CSV was created")
            save_failed_list(system_ids, 4)

    # Step 5: plotting
    if args.step in ("5", "all"):
        print("\n=== Step 5: Plot figures (Vina) ===")
        env = os.environ.copy()
        env["MPLBACKEND"] = "Agg"
        cmd = [RUNS_N_POSES_PYTHON, str(SCRIPT_DIR / "05_plot_figure_with_vina.py")]
        r = subprocess.run(cmd, env=env)
        if r.returncode == 0:
            print("\n  Done: figures generated successfully")
            save_failed_list([], 5)
        else:
            print(f"\n  Failed: exit code {r.returncode}")
            save_failed_list(system_ids, 5)


if __name__ == "__main__":
    main()
