#!/usr/bin/env python
"""
Step 0: End-to-end runner for the docking benchmark on the
runs-n-poses single-ligand dataset.

This script wraps the per-method step 1-4 scripts and is intended to be
launched once in the background; the entire benchmark runs in stages:

  * Ligand preparation with Meeko (writes
    ``meeko_ligands_pdbqt`` and ``meeko_ligands_pdbt``)
  * Docking for each program against the symmetry-corrected inputs
  * Per-method analysis that produces ``predictions/<method>.csv`` and
    ``posebusters_results/<method>.csv``

Existing per-method scripts are reused (not re-implemented) so the
output of this runner is bit-identical to invoking them by hand.

Usage:
    python run_full_benchmark.py --steps ligand,dock,analyze \\
        --threads 8 --max-workers 4 --log-dir logs/benchmark
"""

import argparse
import os
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
DATASETS = Path("/home/rquiroga/Datasets/runs-n-poses-datasets")
PREDICTIONS = DATASETS / "predictions"
POSEBUSTERS = DATASETS / "posebusters_results"
ANALYSIS = Path("/home/rquiroga/github/runs-n-poses/examples/analysis")

# Existing prepared inputs
SYM_REC_PDBQT = DATASETS / "symmetry_receptors_pdbqt"
SYM_REC_PDBT = DATASETS / "symmetry_receptors_pdbt"
SYM_LIG_PDBQT = DATASETS / "symmetry_ligands_pdbqt"
SYM_LIG_PDBT = DATASETS / "symmetry_ligands_pdbt"

# New Meeko-prepared receptors and ligands
MEEKO_REC_PDBQT = DATASETS / "meeko_receptors_pdbqt"
MEEKO_REC_PDBT = DATASETS / "meeko_receptors_pdbt"
MEEKO_LIG_PDBQT = DATASETS / "meeko_ligands_pdbqt"
MEEKO_LIG_PDBT = DATASETS / "meeko_ligands_pdbt"

# Output dirs (per-method docking results)
OUT_VINA_8 = DATASETS / "autodock_vina_8"
OUT_VINA_32 = DATASETS / "autodock_vina_32"
OUT_QVINA_W = DATASETS / "qvina_w"
OUT_QUICKVINA2 = DATASETS / "quickvina2"
OUT_VINARDO = DATASETS / "vinardo_outputs"
OUT_AUTODOCK_GPU = DATASETS / "autodock_gpu"
OUT_RDOCK = DATASETS / "rdock"
OUT_GNINA = DATASETS / "gnina"
OUT_VINA_8_MEEKO = DATASETS / "vina_8_meeko"
OUT_VINA_32_MEEKO = DATASETS / "vina_32_meeko"
OUT_VINARDO_MEEKO = DATASETS / "vinardock_meeko"

PROJECT_PY = "/home/rquiroga/anaconda3/envs/runs_n_poses/bin/python"
SYSTEM_LIST = SCRIPT_DIR / "single_ligand_systems_symmetry.csv"
SINGLE_LIGAND_LIST = SCRIPT_DIR / "single_ligand_systems_symmetry.txt"
SYSTEM_LIST_TXT = SCRIPT_DIR / "single_ligand_systems_symmetry.txt"


def _now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def run_step(cmd: list, log_path: Path, env=None) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    cmd_str = " ".join(shlex.quote(c) for c in cmd)
    print(f"[{_now()}] $ {cmd_str}")
    with open(log_path, "a") as fh:
        fh.write(f"\n[{_now()}] $ {cmd_str}\n")
    with open(log_path, "a") as fh:
        proc = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env)
        rc = proc.wait()
    return rc


def step_ligand_meeko(log_path: Path) -> int:
    # Use the runs_n_poses conda python (has tqdm + meeko available)
    RUNS_N_POSES_PY = "/home/rquiroga/anaconda3/envs/runs_n_poses/bin/python"
    cmd = [
        RUNS_N_POSES_PY, str(SCRIPT_DIR / "prepare_ligand_meeko.py"),
        "--ground-truth-dir", str(DATASETS / "ground_truth"),
        "--output-pdbt", str(MEEKO_LIG_PDBT),
        "--output-pdbqt", str(MEEKO_LIG_PDBQT),
        "--system-info", str(SYSTEM_LIST),
    ]
    return run_step(cmd, log_path)


def step_receptor_meeko(log_path: Path) -> int:
    # Use the runs_n_poses conda python (has tqdm + meeko available)
    RUNS_N_POSES_PY = "/home/rquiroga/anaconda3/envs/runs_n_poses/bin/python"
    cmd = [
        RUNS_N_POSES_PY, str(SCRIPT_DIR / "prepare_receptor_meeko.py"),
        "--input-dir", str(DATASETS / "symmetry_corrected"),
        "--output-pdbqt", str(MEEKO_REC_PDBQT),
        "--output-pdbt", str(MEEKO_REC_PDBT),
        "--system-list", str(SYSTEM_LIST),
    ]
    return run_step(cmd, log_path)


def step_dock_vina(threads: int, log_path: Path, resume: bool) -> int:
    """Run Vina (exh=8 + exh=32) using the obabel-prepared inputs."""
    rc = 0
    # Use the filtered annotations (single-ligand only)
    filtered_annot = SCRIPT_DIR / "annotations_single_ligand_symmetry.csv"
    for exh, out_dir in [(8, OUT_VINA_8), (32, OUT_VINA_32)]:
        cmd = [
            PROJECT_PY, str(SCRIPT_DIR / "03_run_vina.py"),
            "--receptor-dir", str(SYM_REC_PDBQT),
            "--ligand-dir", str(SYM_LIG_PDBQT),
            "--output-dir", str(out_dir),
            "--exhaustiveness", str(exh),
            "--threads", str(threads),
            "--executable", "vina",
            "--annotations", str(filtered_annot),
        ]
        if resume:
            cmd.append("--resume")
        rc |= run_step(cmd, log_path)
    return rc


def step_dock_qvina_w(threads: int, log_path: Path, resume: bool) -> int:
    cmd = [
        PROJECT_PY, str(SCRIPT_DIR / "run_qvina_w.py"),
        "--receptor-dir", str(SYM_REC_PDBQT),
        "--ligand-dir", str(SYM_LIG_PDBQT),
        "--output-dir", str(OUT_QVINA_W),
        "--system-list", str(SYSTEM_LIST_TXT),
        "--cpu", str(threads),
    ]
    if resume:
        cmd.append("--resume")
    return run_step(cmd, log_path)


def step_dock_quickvina2(threads: int, log_path: Path, resume: bool) -> int:
    cmd = [
        PROJECT_PY, str(SCRIPT_DIR / "run_quickvina2.py"),
        "--receptor-dir", str(SYM_REC_PDBQT),
        "--ligand-dir", str(SYM_LIG_PDBQT),
        "--output-dir", str(OUT_QUICKVINA2),
        "--system-list", str(SYSTEM_LIST_TXT),
        "--cpu", str(threads),
    ]
    if resume:
        cmd.append("--resume")
    return run_step(cmd, log_path)


def step_dock_vinardo(threads: int, log_path: Path, resume: bool) -> int:
    cmd = [
        PROJECT_PY, str(SCRIPT_DIR / "run_vinardo_symmetry.py"),
        "--receptor-dir", str(SYM_REC_PDBT),
        "--ligand-dir", str(SYM_LIG_PDBT),
        "--output-dir", str(OUT_VINARDO),
        "--system-list", str(SYSTEM_LIST_TXT),
        "--threads", str(threads),
    ]
    if resume:
        cmd.append("--resume")
    return run_step(cmd, log_path)


def step_dock_autodock_gpu(threads: int, log_path: Path, resume: bool) -> int:
    """AutoDock-GPU: generates PDBQT+GPF from symmetry PDB, docks Meeko ligand."""
    cmd = [
        PROJECT_PY, str(SCRIPT_DIR / "run_autodock_gpu.py"),
        "--sym-receptor-dir", str(DATASETS / "symmetry_corrected"),
        "--ligand-dir", str(MEEKO_LIG_PDBQT),
        "--output-dir", str(OUT_AUTODOCK_GPU),
        "--system-list", str(SYSTEM_LIST_TXT),
        "--nrun", "20",
    ]
    if resume:
        cmd.append("--resume")
    return run_step(cmd, log_path)


def step_dock_vina_meeko_8(threads: int, log_path: Path, resume: bool) -> int:
    """Vina (exh=8) with Meeko-prepared receptor + Meeko-prepared ligand."""
    filtered_annot = SCRIPT_DIR / "annotations_single_ligand_symmetry.csv"
    cmd = [
        PROJECT_PY, str(SCRIPT_DIR / "03_run_vina.py"),
        "--receptor-dir", str(MEEKO_REC_PDBQT),
        "--ligand-dir", str(MEEKO_LIG_PDBQT),
        "--output-dir", str(OUT_VINA_8_MEEKO),
        "--exhaustiveness", "8",
        "--threads", str(threads),
        "--executable", "vina",
        "--annotations", str(filtered_annot),
    ]
    if resume:
        cmd.append("--resume")
    return run_step(cmd, log_path)


def step_dock_vina_meeko_32(threads: int, log_path: Path, resume: bool) -> int:
    """Vina (exh=32) with Meeko-prepared receptor + Meeko-prepared ligand."""
    filtered_annot = SCRIPT_DIR / "annotations_single_ligand_symmetry.csv"
    cmd = [
        PROJECT_PY, str(SCRIPT_DIR / "03_run_vina.py"),
        "--receptor-dir", str(MEEKO_REC_PDBQT),
        "--ligand-dir", str(MEEKO_LIG_PDBQT),
        "--output-dir", str(OUT_VINA_32_MEEKO),
        "--exhaustiveness", "32",
        "--threads", str(threads),
        "--executable", "vina",
        "--annotations", str(filtered_annot),
    ]
    if resume:
        cmd.append("--resume")
    return run_step(cmd, log_path)


def step_dock_gnina(threads: int, log_path: Path, resume: bool) -> int:
    cmd = [
        PROJECT_PY, str(SCRIPT_DIR / "run_gnina.py"),
        "--output-dir", str(OUT_GNINA),
        "--system-list", str(SYSTEM_LIST),
    ]
    if resume:
        cmd.append("--resume")
    return run_step(cmd, log_path)


def step_dock_rdock(threads: int, log_path: Path, resume: bool) -> int:
    cmd = [
        PROJECT_PY, str(SCRIPT_DIR / "run_rdock.py"),
        "--output-dir", str(OUT_RDOCK),
        "--system-list", str(SYSTEM_LIST),
    ]
    if resume:
        cmd.append("--resume")
    return run_step(cmd, log_path)


def step_dock_vinardo_meeko(threads: int, log_path: Path, resume: bool) -> int:
    """vinardock with Meeko-prepared receptor (PDBT) + Meeko-prepared ligand (PDBT)."""
    cmd = [
        PROJECT_PY, str(SCRIPT_DIR / "run_vinardo_symmetry.py"),
        "--receptor-dir", str(MEEKO_REC_PDBT),
        "--ligand-dir", str(MEEKO_LIG_PDBT),
        "--output-dir", str(OUT_VINARDO_MEEKO),
        "--system-list", str(SYSTEM_LIST_TXT),
        "--threads", str(threads),
    ]
    if resume:
        cmd.append("--resume")
    return run_step(cmd, log_path)


# ---------------------------------------------------------------------------
# Analysis wrappers
# ---------------------------------------------------------------------------
def _analyze(family: str, method: str, output_dir: Path, log_path: Path,
             analysis_dir: Path, predictions_csv: Path, posebusters_csv: Path,
             system_list: Path, resume: bool) -> int:
    cmd = [
        PROJECT_PY, str(SCRIPT_DIR / "04_analyze_docking.py"),
        "--family", family,
        "--method", method,
        "--output-dir", str(output_dir),
        "--analysis-dir", str(analysis_dir),
        "--predictions-csv", str(predictions_csv),
        "--posebusters-csv", str(posebusters_csv),
        "--system-list", str(system_list),
    ]
    if resume:
        cmd.append("--resume")
    return run_step(cmd, log_path)


def step_analyze(log_path: Path, resume: bool) -> int:
    rc = 0
    plan = [
        # (family, method, output_dir, analysis_dir, predictions_csv, posebusters_csv, system_list)
        ("vina", "autodock_vina_8", OUT_VINA_8, ANALYSIS / "vina_8",
         PREDICTIONS / "autodock_vina_8.csv", POSEBUSTERS / "autodock_vina_8.csv", SYSTEM_LIST_TXT),
        ("vina", "autodock_vina_32", OUT_VINA_32, ANALYSIS / "vina_32",
         PREDICTIONS / "autodock_vina_32.csv", POSEBUSTERS / "autodock_vina_32.csv", SYSTEM_LIST_TXT),
        ("vina", "qvina_w", OUT_QVINA_W, ANALYSIS / "qvina_w",
         PREDICTIONS / "qvina_w.csv", POSEBUSTERS / "qvina_w.csv", SYSTEM_LIST_TXT),
        ("vina", "quickvina2", OUT_QUICKVINA2, ANALYSIS / "quickvina2",
         PREDICTIONS / "quickvina2.csv", POSEBUSTERS / "quickvina2.csv", SYSTEM_LIST_TXT),
        ("vina", "vina_8_meeko", OUT_VINA_8_MEEKO, ANALYSIS / "vina_8_meeko",
         PREDICTIONS / "vina_8_meeko.csv", POSEBUSTERS / "vina_8_meeko.csv", SYSTEM_LIST_TXT),
        ("vina", "vina_32_meeko", OUT_VINA_32_MEEKO, ANALYSIS / "vina_32_meeko",
         PREDICTIONS / "vina_32_meeko.csv", POSEBUSTERS / "vina_32_meeko.csv", SYSTEM_LIST_TXT),
        ("vinardo", "vinardock_2vinardo", OUT_VINARDO, ANALYSIS / "vinardock_2vinardo",
         PREDICTIONS / "vinardock_2vinardo.csv", POSEBUSTERS / "vinardock_2vinardo.csv", SYSTEM_LIST_TXT),
        ("vinardo", "vinardock_meeko", OUT_VINARDO_MEEKO, ANALYSIS / "vinardock_meeko",
         PREDICTIONS / "vinardock_meeko.csv", POSEBUSTERS / "vinardock_meeko.csv", SYSTEM_LIST_TXT),
        ("vina", "gnina", OUT_GNINA, ANALYSIS / "gnina",
         PREDICTIONS / "gnina.csv", POSEBUSTERS / "gnina.csv", SYSTEM_LIST_TXT),
        ("vina", "rdock", OUT_RDOCK, ANALYSIS / "rdock",
         PREDICTIONS / "rdock.csv", POSEBUSTERS / "rdock.csv", SYSTEM_LIST_TXT),
        ("adgpu", "autodock_gpu", OUT_AUTODOCK_GPU, ANALYSIS / "autodock_gpu",
         PREDICTIONS / "autodock_gpu.csv", POSEBUSTERS / "autodock_gpu.csv", SYSTEM_LIST_TXT),
        # autodock_gpu_mostpop is analyzed inline inside the adgpu analysis (same output dir)
    ]
    for family, method, out_dir, analysis_dir, pred_csv, pb_csv, sys_list in plan:
        rc |= _analyze(family, method, out_dir, log_path,
                       analysis_dir, pred_csv, pb_csv, sys_list, resume)
    return rc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", default="receptor,ligand,dock,analyze",
                        help="comma-separated stages to run")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--max-workers", type=int, default=1,
                        help="placeholder for future parallelization; methods run sequentially")
    parser.add_argument("--log-dir", default=str(SCRIPT_DIR / "logs/benchmark"),
                        help="directory for combined run logs")
    parser.add_argument("--resume", action="store_true",
                        help="skip already-completed per-method outputs")
    args = parser.parse_args()

    log_dir = Path(args.log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    main_log = log_dir / "benchmark.log"

    print(f"Logging to {main_log}")
    with open(main_log, "a") as fh:
        fh.write(f"\n\n==== run_full_benchmark @ {_now()} ====\n")

    steps = [s.strip() for s in args.steps.split(",") if s.strip()]

    if "receptor" in steps:
        print("\n=== Stage: Meeko receptor prep ===")
        rc = step_receptor_meeko(main_log)
        if rc != 0:
            print(f"  Meeko receptor prep exited with {rc}")
            sys.exit(rc)

    if "ligand" in steps:
        print("\n=== Stage: Meeko ligand prep ===")
        rc = step_ligand_meeko(main_log)
        if rc != 0:
            print(f"  Meeko ligand prep exited with {rc}")
            sys.exit(rc)

    if "dock" in steps or "analyze" in steps:
        print("\n=== Stage: Dock + Analyze ===")

        # Simple: for each method, dock (skipped if already done by --resume), then analyze.
        methods = [
            ("vina (8,32)",    step_dock_vina,        [("vina", "autodock_vina_8",   OUT_VINA_8,    ANALYSIS / "vina_8"),
                                                       ("vina", "autodock_vina_32",  OUT_VINA_32,   ANALYSIS / "vina_32")]),
            ("qvina-w",        step_dock_qvina_w,     [("vina", "qvina_w",           OUT_QVINA_W,   ANALYSIS / "qvina_w")]),
            ("quickvina2",     step_dock_quickvina2,  [("vina", "quickvina2",        OUT_QUICKVINA2,ANALYSIS / "quickvina2")]),
            ("autodock-gpu",   step_dock_autodock_gpu,[("adgpu","autodock_gpu",      OUT_AUTODOCK_GPU,ANALYSIS / "autodock_gpu")]),
            ("vina-8-meeko",   step_dock_vina_meeko_8,[("vina", "vina_8_meeko",      OUT_VINA_8_MEEKO,ANALYSIS / "vina_8_meeko")]),
            ("vinardock",      step_dock_vinardo,     [("vinardo","vinardock_2vinardo",OUT_VINARDO,ANALYSIS / "vinardock_2vinardo")]),
            ("vinardock-meeko",step_dock_vinardo_meeko,[("vinardo","vinardock_meeko",OUT_VINARDO_MEEKO,ANALYSIS / "vinardock_meeko")]),
            ("rdock",          step_dock_rdock,       [("vina", "rdock",             OUT_RDOCK,     ANALYSIS / "rdock")]),
            ("gnina",          step_dock_gnina,       [("vina", "gnina",             OUT_GNINA,     ANALYSIS / "gnina")]),
            ("vina-32-meeko",  step_dock_vina_meeko_32,[("vina","vina_32_meeko",     OUT_VINA_32_MEEKO,ANALYSIS / "vina_32_meeko")]),
        ]

        for label, dock_fn, analyze_list in methods:
            if "dock" in steps and dock_fn is not None:
                print(f"\n--- Dock: {label} ---")
                dock_fn(args.threads, main_log, args.resume)
            if "analyze" in steps:
                for family, method, out_dir, analysis_dir in analyze_list:
                    print(f"\n--- Analyze: {method} ---")
                    _analyze(family, method, out_dir, main_log,
                             analysis_dir, PREDICTIONS / f"{method}.csv",
                             POSEBUSTERS / f"{method}.csv", SYSTEM_LIST_TXT, args.resume)

    print(f"\nAll done. Combined log: {main_log}")


if __name__ == "__main__":
    main()
