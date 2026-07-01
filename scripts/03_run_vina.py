#!/usr/bin/env python
"""
Step 3 (VINA): Run AutoDock Vina for each complex.

This script computes a simple autobox from receptor coordinates and runs `vina` for
each receptor-ligand pair. Output layout mirrors vinardock outputs: a per-system
folder with subfolder `<system>_<ligand>` containing `out.pdbqt` and `log.txt`.

Usage:
    python 03_run_vina.py --receptor-dir /path/to/receptors \
                         --ligand-dir /path/to/ligands \
                         --output-dir /path/to/vina_outputs
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path
from tqdm import tqdm
import math
import time
import json


def parse_pdb_coordinates(pdbfile: str):
    """
    Parse atom coordinates from a PDB/PDBQT/PDBT file by reading ATOM/HETATM lines.
    Returns list of (x,y,z) tuples.
    """
    coords = []
    try:
        with open(pdbfile, 'r') as fh:
            for line in fh:
                if line.startswith(('ATOM', 'HETATM')) and len(line) >= 54:
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


def compute_autobox_from_receptor(receptor_file: str, margin: float = 10.0):
    coords = parse_pdb_coordinates(receptor_file)
    if not coords:
        return None
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    zmin, zmax = min(zs), max(zs)
    cx = (xmin + xmax) / 2.0
    cy = (ymin + ymax) / 2.0
    cz = (zmin + zmax) / 2.0
    sx = max((xmax - xmin) + margin, 5.0)
    sy = max((ymax - ymin) + margin, 5.0)
    sz = max((zmax - zmin) + margin, 5.0)
    return cx, cy, cz, sx, sy, sz


def compute_autobox_from_ligand(ligand_file: str, margin: float = 22.0):
    """Compute box from ligand coordinates (much smaller, avoids huge boxes)."""
    coords = parse_pdb_coordinates(ligand_file)
    if not coords:
        return None
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    cx = (min(xs) + max(xs)) / 2.0
    cy = (min(ys) + max(ys)) / 2.0
    cz = (min(zs) + max(zs)) / 2.0
    sx = max((max(xs) - min(xs)) + margin, 22.0)
    sy = max((max(ys) - min(ys)) + margin, 22.0)
    sz = max((max(zs) - min(zs)) + margin, 22.0)
    return cx, cy, cz, sx, sy, sz


def run_vina(receptor_file: str, ligand_file: str, output_dir: str,
             system_id: str, ligand_chain: str, executable: str = "vina",
             threads: int = 4, exhaustiveness: int = 8, margin: float = 10.0,
             box_from_ligand: bool = False,
             box_cx: float = None, box_cy: float = None, box_cz: float = None,
             box_sx: float = None, box_sy: float = None, box_sz: float = None) -> bool:
    system_output_dir = Path(output_dir) / system_id
    system_output_dir.mkdir(parents=True, exist_ok=True)

    out_folder = system_output_dir / f"{system_id}_{ligand_chain}"
    out_folder.mkdir(exist_ok=True)

    # output filenames
    out_pdbqt = out_folder / "out.pdbqt"
    log_file = out_folder / "log.txt"

    # If receptor PDBQT contains non-standard tags (e.g., ROOT) that Vina rejects,
    # attempt a round-trip conversion via PDB to clean the file.
    try:
        with open(receptor_file, 'r') as fh:
            content = fh.read()
        if 'ROOT' in content:
            cleaned = str(Path(receptor_file).with_name(Path(receptor_file).stem + '_clean.pdbqt'))
            # Round-trip: pdbqt -> pdb -> pdbqt
            try:
                subprocess.run(["obabel-25-07", receptor_file, "-O", cleaned.replace('_clean.pdbqt', '_roundtrip.pdb')], capture_output=True, text=True, timeout=60)
                subprocess.run(["obabel-25-07", cleaned.replace('_clean.pdbqt', '_roundtrip.pdb'), "-O", cleaned, "-h", "-xc", "-xr"], capture_output=True, text=True, timeout=60)
                # If cleaned file created, switch to using it
                if os.path.exists(cleaned) and os.path.getsize(cleaned) > 0:
                    receptor_file = cleaned
            except Exception:
                pass
    except Exception:
        pass

    if box_cx is not None:
        cx, cy, cz, sx, sy, sz = box_cx, box_cy, box_cz, box_sx, box_sy, box_sz
    elif box_from_ligand:
        box = compute_autobox_from_ligand(ligand_file, margin=22.0)
        if box is None:
            print(f"Warning: failed to compute box for {system_id}")
            return False
        cx, cy, cz, sx, sy, sz = box
    else:
        box = compute_autobox_from_receptor(receptor_file, margin=margin)
        if box is None:
            print(f"Warning: failed to compute box for {system_id}")
            return False
        cx, cy, cz, sx, sy, sz = box

    cmd = [
        executable,
        "--receptor",
        receptor_file,
        "--ligand",
        ligand_file,
        "--out",
        str(out_pdbqt),
        "--cpu",
        str(threads),
        "--exhaustiveness",
        str(exhaustiveness),
        "--center_x",
        str(cx),
        "--center_y",
        str(cy),
        "--center_z",
        str(cz),
        "--size_x",
        str(sx),
        "--size_y",
        str(sy),
        "--size_z",
        str(sz),
    ]

    try:
        start_time = time.time()
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=900,
            cwd=str(system_output_dir),
        )
        end_time = time.time()
        elapsed = end_time - start_time

        # Save stdout/stderr to log_file for diagnostics
        try:
            with open(log_file, 'w') as fh:
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
                "method": "vina",
                "exhaustiveness": exhaustiveness,
                "returncode": result.returncode,
                "runtime_seconds": float(elapsed),
            }
            runtime_file = out_folder / "runtime.json"
            with open(runtime_file, 'w') as fh:
                json.dump(run_meta, fh)
        except Exception:
            pass

        if result.returncode != 0:
            print(f"Warning: vina failed for {system_id} {ligand_chain}:")
            if result.stderr:
                print(f"  stderr: {result.stderr.splitlines()[0]}")
            return False

        # Check output existence
        if out_pdbqt.exists() and out_pdbqt.stat().st_size > 0:
            return True
        print(f"Warning: No output files created for {system_id} {ligand_chain}")
        return False

    except subprocess.TimeoutExpired:
        print(f"Warning: vina timed out for {system_id} {ligand_chain}")
        return False
    except Exception as e:
        print(f"Error running vina for {system_id} {ligand_chain}: {e}")
        return False


def process_system(system_id: str, receptor_dir: Path, ligand_dir: Path,
                   output_dir: Path, executable: str, threads: int = 4,
                   exhaustiveness: int = 8, resume: bool = False,
                   box_from_ligand: bool = False,
                   box_cx: float = None, box_cy: float = None, box_cz: float = None,
                   box_sx: float = None, box_sy: float = None, box_sz: float = None) -> bool:
    receptor_file = receptor_dir / system_id / f"{system_id}_receptor.pdbqt"
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

    ligand_files = list(ligand_dir_system.glob("*.pdbqt"))
    if not ligand_files:
        print(f"Warning: No ligand files found for {system_id}")
        return False

    success = False
    for ligand_file in ligand_files:
        ligand_chain = ligand_file.stem
        print(f"Running vina for {system_id}, ligand {ligand_chain}...")

        out_folder = Path(output_dir) / system_id / f"{system_id}_{ligand_chain}"
        out_pdbqt = out_folder / "out.pdbqt"
        log_file = out_folder / "log.txt"

        if resume and out_pdbqt.exists() and out_pdbqt.stat().st_size > 0:
            print(f"  Skipping (exists) {system_id} {ligand_chain}")
            success = True
            continue
        if resume and log_file.exists() and log_file.stat().st_size > 50:
            print(f"  Skipping (log exists) {system_id} {ligand_chain}")
            success = True
            continue

        if run_vina(
            str(receptor_file),
            str(ligand_file),
            str(output_dir),
            system_id,
            ligand_chain,
            executable,
            threads,
            exhaustiveness,
            box_from_ligand=box_from_ligand,
            box_cx=box_cx, box_cy=box_cy, box_cz=box_cz,
            box_sx=box_sx, box_sy=box_sy, box_sz=box_sz,
        ):
            success = True
        else:
            print(f"  Failed for ligand {ligand_chain}")

    return success


def main():
    parser = argparse.ArgumentParser(description="Run AutoDock Vina for each complex")
    parser.add_argument(
        "--receptor-dir",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/vina_inputs/receptors",
    )
    parser.add_argument(
        "--ligand-dir",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/vina_inputs/ligands",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/vina_outputs",
    )
    parser.add_argument(
        "--executable",
        type=str,
        default="vina",
    )
    parser.add_argument(
        "--exhaustiveness",
        type=int,
        default=8,
        help="Vina exhaustiveness value (higher -> more thorough search)",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Skip ligand runs when output exists (out.pdbqt non-empty or log.txt present)",
    )
    parser.add_argument("--system-id", type=str, default=None)
    parser.add_argument("--annotations", type=str, default="/home/rquiroga/Datasets/runs-n-poses-datasets/annotations.csv")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--center-x", type=float, default=None, help="Box center X (overrides autobox)")
    parser.add_argument("--center-y", type=float, default=None)
    parser.add_argument("--center-z", type=float, default=None)
    parser.add_argument("--size-x", type=float, default=None, help="Box size X (overrides autobox)")
    parser.add_argument("--size-y", type=float, default=None)
    parser.add_argument("--size-z", type=float, default=None)

    args = parser.parse_args()

    receptor_dir = Path(args.receptor_dir)
    ligand_dir = Path(args.ligand_dir)
    # allow comma-separated list of exhaustiveness values
    exh_values = [int(x) for x in str(args.exhaustiveness).split(",") if str(x).strip()]

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if os.path.exists(args.annotations):
        import pandas as pd
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

    # Load pre-computed boxes
    boxes_path = Path(__file__).parent / "boxes.json"
    if boxes_path.exists():
        import json
        boxes = json.loads(boxes_path.read_text())
    else:
        boxes = {}

    success_count = 0
    fail_count = 0

    for exh in exh_values:
        if len(exh_values) == 1:
            out_dir_ex = output_dir
        else:
            out_dir_ex = output_dir.parent / f"{output_dir.name}_{exh}"
        out_dir_ex.mkdir(parents=True, exist_ok=True)

        print(f"Running Vina (exhaustiveness={exh}) -> {out_dir_ex}")

        for system_id in tqdm(system_ids, desc=f"Running Vina (ex={exh})"):
            b = boxes.get(system_id, {})
            if process_system(
                system_id,
                receptor_dir,
                ligand_dir,
                out_dir_ex,
                args.executable,
                args.threads,
                exh,
                resume=args.resume,
                box_cx=b.get("center_x"), box_cy=b.get("center_y"), box_cz=b.get("center_z"),
                box_sx=b.get("size_x"), box_sy=b.get("size_y"), box_sz=b.get("size_z"),
            ):
                success_count += 1
            else:
                fail_count += 1

    print(f"\nDone! Success: {success_count}, Failed: {fail_count}")
    print(f"Output directory (base): {output_dir}")


if __name__ == "__main__":
    main()
