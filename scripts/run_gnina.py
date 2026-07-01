#!/usr/bin/env python
"""
Run GNINA on symmetry-corrected receptors.

GNINA is a deep learning-augmented docking tool based on AutoDock Vina.
It uses CNN scoring to rank poses. Requires LD_PRELOAD for the custom
OpenBabel fork that provides SetSegName.

Usage:
    python run_gnina.py \
        --output-dir runs-n-poses-datasets/gnina \
        --system-list scripts/single_ligand_systems_symmetry.csv
"""

import argparse
import csv
import json
import os
import subprocess
import time
from pathlib import Path

GNINA = "/home/rquiroga/.local/bin/gnina"
OBABEL_LIB = "/home/rquiroga/Downloads/openbabel-julio25/openbabel-master/build/lib/libopenbabel.so.7"
SYM_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/symmetry_corrected")
GT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth")


def parse_coordinates(filepath):
    """Parse 3D coordinates from PDB, SDF or PDBQT files."""
    coords = []
    with open(filepath) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                try:
                    coords.append((
                        float(line[30:38].strip()),
                        float(line[38:46].strip()),
                        float(line[46:54].strip()),
                    ))
                except Exception:
                    continue
    return coords


def run_gnina(sys_id: str, chain: str, output_dir: Path, resume: bool) -> tuple[str, float | None]:
    out_dir = output_dir / sys_id / f"{sys_id}_{chain}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_sdf = out_dir / "out.sdf"
    log_file = out_dir / "log.txt"
    runtime_file = out_dir / "runtime.json"

    if resume and out_sdf.exists() and out_sdf.stat().st_size > 0:
        return "skip", None

    sym_pdb = SYM_DIR / sys_id / f"{sys_id}_receptor_symm.pdb"
    if not sym_pdb.exists():
        return "no_receptor", None

    lig_sdf = GT_DIR / sys_id / "ligand_files" / f"{chain}.sdf"
    if not lig_sdf.exists():
        return "no_ligand", None

    # Autobox from ligand coordinates (convert SDF→PDB for parsing)
    lig_pdb = out_dir / "ligand.pdb"
    subprocess.run(["obabel-25-07", str(lig_sdf), "-O", str(lig_pdb), "-h"],
                   capture_output=True, timeout=30)
    coords = parse_coordinates(str(lig_pdb)) if lig_pdb.exists() else []
    if not coords:
        return "no_coords", None

    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    cx = (min(xs) + max(xs)) / 2.0
    cy = (min(ys) + max(ys)) / 2.0
    cz = (min(zs) + max(zs)) / 2.0
    sx = max(max(xs) - min(xs) + 15.0, 20.0)
    sy = max(max(ys) - min(ys) + 15.0, 20.0)
    sz = max(max(zs) - min(zs) + 15.0, 20.0)

    env = os.environ.copy()
    env["LD_PRELOAD"] = OBABEL_LIB + ":" + env.get("LD_PRELOAD", "")

    cmd = [
        GNINA,
        "-r", str(sym_pdb),
        "-l", str(lig_sdf),
        "--autobox_ligand", str(lig_sdf),
        "--center_x", f"{cx:.3f}",
        "--center_y", f"{cy:.3f}",
        "--center_z", f"{cz:.3f}",
        "--size_x", f"{sx:.3f}",
        "--size_y", f"{sy:.3f}",
        "--size_z", f"{sz:.3f}",
        "--num_modes", "10",
        "--cpu", "8",
        "-o", str(out_sdf),
    ]

    try:
        start = time.time()
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600, env=env)
        elapsed = time.time() - start

        with open(log_file, "w") as f:
            if r.stdout:
                f.write(r.stdout)
            if r.stderr:
                f.write("\nSTDERR:\n" + r.stderr)

        if r.returncode == 0 and out_sdf.exists() and out_sdf.stat().st_size > 0:
            with open(runtime_file, "w") as f:
                json.dump({"method": "gnina", "runtime_seconds": elapsed}, f)
            return "ok", elapsed

        return "dock_fail", None
    except subprocess.TimeoutExpired:
        return "timeout", None
    except Exception as e:
        return f"error: {e}", None


def main():
    parser = argparse.ArgumentParser(description="Run GNINA on symmetry-corrected receptors")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--system-list", required=True)
    parser.add_argument("--system-id", type=str, default=None)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(args.system_list) as f:
        reader = csv.DictReader(f)
        systems = [(r["system_id"], r["proper_ligand_chain"]) for r in reader]
    if args.system_id:
        systems = [s for s in systems if s[0] == args.system_id]

    ok, skip, fail = 0, 0, 0
    for sys_id, chain in systems:
        status, rt = run_gnina(sys_id, chain, output_dir, args.resume)
        if status == "ok":
            print(f"  OK {sys_id} ({rt:.1f}s)")
            ok += 1
        elif status == "skip":
            skip += 1
        else:
            print(f"  {status}: {sys_id}")
            fail += 1

    print(f"GNINA: OK={ok}, Skip={skip}, Fail={fail}")


if __name__ == "__main__":
    main()
