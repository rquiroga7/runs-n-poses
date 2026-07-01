#!/usr/bin/env python
"""
Run LeDock on symmetry-corrected receptors.

LeDock uses a simple config file with receptor PDB, ligand MOL2/SDF,
and box coordinates derived from the input ligand.

Usage:
    python run_ledock.py \
        --sym-dir runs-n-poses-datasets/symmetry_corrected \
        --ground-truth runs-n-poses-datasets/ground_truth \
        --output-dir runs-n-poses-datasets/ledock \
        --system-list scripts/single_ligand_systems_symmetry.csv
"""

import argparse
import csv
import os
import subprocess
import tempfile
import time
from pathlib import Path

from rdkit import Chem

LEDOCK = "/home/rquiroga/Downloads/ledock_linux_x86"
SYM_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/symmetry_corrected")
GT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth")


def sdf_to_mol2(sdf_path: str, mol2_path: str) -> bool:
    """Convert SDF to MOL2 using obabel."""
    try:
        r = subprocess.run(
            ["obabel-25-07", sdf_path, "-O", mol2_path, "-h"],
            capture_output=True, text=True, timeout=30,
        )
        return r.returncode == 0 and os.path.exists(mol2_path)
    except Exception:
        return False


def write_config(cfg_path: str, receptor_pdb: str, ligand_mol2: str,
                 output_dok: str, cx: float, cy: float, cz: float,
                 rx: float, ry: float, rz: float, n_runs: int = 100) -> None:
    with open(cfg_path, "w") as f:
        f.write(f"Receptor: {receptor_pdb}\n")
        f.write(f"Output: {output_dok}\n")
        f.write(f"Ligand: {ligand_mol2}\n")
        f.write(f"Box center: {cx:.3f} {cy:.3f} {cz:.3f}\n")
        f.write(f"Box radius: {rx:.3f} {ry:.3f} {rz:.3f}\n")
        f.write(f"Number of runs: {n_runs}\n")


def parse_pdb_coordinates(pdbfile: str):
    coords = []
    with open(pdbfile) as fh:
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


def run_ledock(sys_id: str, chain: str, output_dir: Path, resume: bool) -> tuple[str, float | None]:
    out_dir = output_dir / sys_id / f"{sys_id}_{chain}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_dok = out_dir / "out.dok"
    log_file = out_dir / "log.txt"
    runtime_file = out_dir / "runtime.json"

    if resume and out_dok.exists() and out_dok.stat().st_size > 0:
        return "skip", None

    sym_pdb = SYM_DIR / sys_id / f"{sys_id}_receptor_symm.pdb"
    if not sym_pdb.exists():
        return "no_receptor", None

    lig_sdf = GT_DIR / sys_id / "ligand_files" / f"{chain}.sdf"
    if not lig_sdf.exists():
        return "no_ligand", None

    # Convert SDF to MOL2 for LeDock
    lig_mol2 = out_dir / "ligand.mol2"
    if not lig_mol2.exists() or lig_mol2.stat().st_size == 0:
        if not sdf_to_mol2(str(lig_sdf), str(lig_mol2)):
            return "mol2_fail", None

    # Box from ligand coords
    coords = parse_pdb_coordinates(str(sym_pdb))
    if not coords:
        return "no_coords", None
    # Use ligand coordinates for box center
    lig_coords = parse_pdb_coordinates(str(lig_mol2))
    if not lig_coords:
        lig_coords = parse_pdb_coordinates(str(lig_sdf))
    if lig_coords:
        xs = [c[0] for c in lig_coords]
        ys = [c[1] for c in lig_coords]
        zs = [c[2] for c in lig_coords]
        cx = (min(xs) + max(xs)) / 2.0
        cy = (min(ys) + max(ys)) / 2.0
        cz = (min(zs) + max(zs)) / 2.0
        rx = max(max(xs) - min(xs) + 15.0, 20.0) / 2.0
        ry = max(max(ys) - min(ys) + 15.0, 20.0) / 2.0
        rz = max(max(zs) - min(zs) + 15.0, 20.0) / 2.0
    else:
        # Fallback to receptor center
        xs = [c[0] for c in coords]
        ys = [c[1] for c in coords]
        zs = [c[2] for c in coords]
        cx = (min(xs) + max(xs)) / 2.0
        cy = (min(ys) + max(ys)) / 2.0
        cz = (min(zs) + max(zs)) / 2.0
        rx = max(xs) - min(xs) + 10.0
        ry = max(ys) - min(ys) + 10.0
        rz = max(zs) - min(zs) + 10.0

    cfg_file = out_dir / "ledock.cfg"
    write_config(str(cfg_file), str(sym_pdb), str(lig_mol2),
                 str(out_dok), cx, cy, cz, rx, ry, rz)

    try:
        start = time.time()
        result = subprocess.run(
            [LEDOCK, str(cfg_file)],
            capture_output=True, text=True, timeout=600,
            cwd=str(out_dir),
        )
        elapsed = time.time() - start

        with open(log_file, "w") as f:
            if result.stdout:
                f.write(result.stdout)
            if result.stderr:
                f.write("\nSTDERR:\n" + result.stderr)

        if result.returncode == 0 and out_dok.exists() and out_dok.stat().st_size > 0:
            with open(runtime_file, "w") as f:
                import json
                json.dump({"method": "ledock", "runtime_seconds": elapsed}, f)
            return "ok", elapsed

        return "dock_fail", None
    except subprocess.TimeoutExpired:
        return "timeout", None
    except Exception as e:
        return f"error: {e}", None


def main():
    parser = argparse.ArgumentParser(description="Run LeDock on symmetry-corrected receptors")
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
        status, rt = run_ledock(sys_id, chain, output_dir, args.resume)
        if status == "ok":
            print(f"  OK {sys_id} ({rt:.1f}s)")
            ok += 1
        elif status == "skip":
            skip += 1
        else:
            print(f"  {status}: {sys_id}")
            fail += 1

    print(f"LeDock: OK={ok}, Skip={skip}, Fail={fail}")


if __name__ == "__main__":
    main()
