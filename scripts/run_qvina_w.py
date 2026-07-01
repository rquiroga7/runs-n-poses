#!/usr/bin/env python
"""
Step 4c: Run QuickVina-W on symmetry-corrected receptors.

QuickVina-W is a fork of QuickVina2 focused on wider binding sites and
parallelism. It accepts standard Vina PDBQT inputs.

Usage:
    python run_qvina_w.py \
        --receptor-dir runs-n-poses-datasets/symmetry_receptors_pdbqt \
        --ligand-dir runs-n-poses-datasets/symmetry_ligands_pdbqt \
        --output-dir runs-n-poses-datasets/qvina_w_symmetry \
        --system-list scripts/systems_for_symmetry_docking.txt
"""

import argparse
import json
import os
import subprocess
import time
from pathlib import Path

QVINA_W = "qvina-w"


def parse_pdb_coordinates(pdbfile):
    """Parse atom coordinates from a PDB/PDBQT file."""
    coords = []
    try:
        with open(pdbfile) as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        coords.append((x, y, z))
                    except Exception:
                        continue
    except Exception:
        pass
    return coords


def compute_autobox(coords, margin=10.0):
    """Compute box center and dimensions from coordinates + margin."""
    if not coords:
        return None
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    cx = (min(xs) + max(xs)) / 2.0
    cy = (min(ys) + max(ys)) / 2.0
    cz = (min(zs) + max(zs)) / 2.0
    sx = max(max(xs) - min(xs) + margin, 5.0)
    sy = max(max(ys) - min(ys) + margin, 5.0)
    sz = max(max(zs) - min(zs) + margin, 5.0)
    return cx, cy, cz, sx, sy, sz


def run_qvina_w(receptor_file, ligand_file, output_dir, sys_id, lig_chain,
                cpu=8, exhaustiveness=8, margin=10.0):
    """Run QuickVina-W docking for a single receptor-ligand pair."""
    output_dir = Path(output_dir).resolve()
    out_folder = output_dir / sys_id / f"{sys_id}_{lig_chain}"
    out_folder.mkdir(parents=True, exist_ok=True)

    out_pdbqt = out_folder / "out.pdbqt"
    log_file = out_folder / "log.txt"

    box = compute_autobox(parse_pdb_coordinates(ligand_file), margin=22.0)
    if box is None:
        return False
    cx, cy, cz, sx, sy, sz = box

    cmd = [
        QVINA_W,
        "--receptor", str(receptor_file),
        "--ligand", str(ligand_file),
        "--out", str(out_pdbqt),
        "--cpu", str(cpu),
        "--exhaustiveness", str(exhaustiveness),
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
    ]

    try:
        start = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=900,
                                cwd=str(out_folder.parent))
        elapsed = time.time() - start

        with open(log_file, "w") as fh:
            if result.stdout:
                fh.write("STDOUT:\n" + result.stdout)
            if result.stderr:
                fh.write("\nSTDERR:\n" + result.stderr)

        runtime_file = out_folder / "runtime.json"
        with open(runtime_file, "w") as fh:
            json.dump({"method": "qvina-w", "exhaustiveness": exhaustiveness,
                       "returncode": result.returncode,
                       "runtime_seconds": elapsed}, fh)

        if result.returncode != 0:
            return False
        return out_pdbqt.exists() and out_pdbqt.stat().st_size > 0
    except subprocess.TimeoutExpired:
        return False
    except Exception as e:
        print(f"Error: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(description="Run QuickVina-W on symmetry-corrected receptors")
    parser.add_argument("--receptor-dir", required=True)
    parser.add_argument("--ligand-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--system-list", required=True)
    parser.add_argument("--system-id", type=str, default=None)
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--cpu", type=int, default=8)
    parser.add_argument("--margin", type=float, default=10.0)
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
        receptor_file = (receptor_dir / sys_id / f"{sys_id}_receptor.pdbqt").resolve()
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

        ligand_files = sorted(lig_dir.glob("*.pdbqt"))
        if not ligand_files:
            print(f"  No ligand files: {sys_id}")
            fail += 1
            continue

        sys_ok = False
        for lf in ligand_files:
            lig_chain = lf.stem
            lf = lf.resolve()
            out_folder = output_dir / sys_id / f"{sys_id}_{lig_chain}"
            out_pdbqt = out_folder / "out.pdbqt"

            if args.resume and out_pdbqt.exists() and out_pdbqt.stat().st_size > 0:
                skip += 1
                sys_ok = True
                continue

            if run_qvina_w(str(receptor_file), str(lf), str(output_dir),
                           sys_id, lig_chain, cpu=args.cpu,
                           exhaustiveness=args.exhaustiveness, margin=args.margin):
                ok += 1
                sys_ok = True
            else:
                print(f"  FAIL qvina-w: {sys_id} {lig_chain}")
                fail += 1

        if not sys_ok:
            fail += 1

    print(f"QuickVina-W: OK={ok}, Skip={skip}, Fail={fail}")


if __name__ == "__main__":
    main()
