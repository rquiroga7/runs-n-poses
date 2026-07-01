#!/usr/bin/env python
"""
Step 4g: Run AutoDock-GPU (Scripps) on symmetry-corrected receptors.

Pipeline per system:
  1. Use the symmetry-corrected receptor PDB -> mk_prepare_receptor
     (``--read_pdb`` / ``-p -g``) to generate PDBQT + GPF in one call.
  2. The Meeko-prepared ligand PDBQT is used directly.
  3. The box is derived from the ligand coordinates.
  4. autogrid4 -> grid maps (.maps.fld)
  5. AutoDock-GPU -> docking results (.dlg)

Usage:
    python run_autodock_gpu.py \
        --sym-receptor-dir runs-n-poses-datasets/symmetry_corrected \
        --ligand-dir runs-n-poses-datasets/meeko_ligands_pdbqt \
        --output-dir runs-n-poses-datasets/autodock_gpu \
        --system-list scripts/single_ligand_systems_symmetry.csv

Note: systems are processed sequentially so per-system runtime is accurate.
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

ADGPU = "AutoDock-GPU"
AUTOGRID4 = "/home/rquiroga/Downloads/autodocksuite-4.2.6-x86_64Linux2/x86_64Linux2/autogrid4"
MEEKO_REC = "/home/rquiroga/github/runs-n-poses/.venv/bin/python -m meeko.cli.mk_prepare_receptor"


def get_receptor_types(pdbqt_path):
    types = set()
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                parts = line.split()
                if len(parts) >= 12:
                    types.add(parts[-1])
    return sorted(types)


def get_ligand_types(pdbqt_path):
    types = set()
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                parts = line.split()
                if len(parts) >= 12:
                    types.add(parts[-1])
    return sorted(types)


def fix_gpf(gpf_path, rec_types, lig_types):
    """Rewrite receptor_types / ligand_types / map lines in the GPF."""
    with open(gpf_path) as f:
        lines = f.readlines()

    new_lines = []
    seen_rec = False
    for line in lines:
        if line.startswith("receptor_types"):
            new_lines.append("receptor_types " + " ".join(rec_types) + "\n")
            seen_rec = True
        else:
            new_lines.append(line)
    if not seen_rec:
        raise RuntimeError(f"No receptor_types line in {gpf_path}")

    gpf = "".join(new_lines)
    gpf = re.sub(
        r"^ligand_types .*",
        "ligand_types " + " ".join(lig_types),
        gpf,
        flags=re.MULTILINE,
    )
    gpf = re.sub(
        r"^map .*?\.(\w+)\.map",
        lambda m: m.group(0) if m.group(1) in lig_types else "",
        gpf,
        flags=re.MULTILINE,
    )
    gpf = re.sub(r"\n{2,}", "\n", gpf).strip() + "\n"

    with open(gpf_path, "w") as f:
        f.write(gpf)


def run_autogrid4(gpf_path, glg_path):
    result = subprocess.run(
        [AUTOGRID4, "-p", str(gpf_path), "-l", str(glg_path)],
        capture_output=True, text=True, timeout=300,
        cwd=str(gpf_path.parent),
    )
    return result.returncode == 0


def _extract_mostpop_pose(dlg_path: Path, out_pdbqt: Path) -> bool:
    """Extract the best-energy pose from the most-populated cluster in a DLG."""
    try:
        with open(dlg_path) as f:
            text = f.read()
        import re
        hist_match = re.search(r'CLUSTERING HISTOGRAM.*?\n(.*?)(?=\n\s+RMSD TABLE)', text, re.DOTALL)
        if not hist_match:
            return False
        hist = hist_match.group(1)
        clusters = []
        for line in hist.strip().split('\n'):
            line = line.strip()
            if not line or '|' not in line or line.startswith('Clus') or '____' in line:
                continue
            parts = line.split('|')
            if len(parts) < 6:
                continue
            try:
                rank = int(parts[0].strip())
                energy = float(parts[1].strip())
                run = int(parts[2].strip())
                mean_energy = float(parts[3].strip())
                num_in_clus = int(parts[4].strip())
                clusters.append((rank, energy, run, mean_energy, num_in_clus))
            except (ValueError, IndexError):
                continue
        if not clusters:
            return False
        # Most populated cluster (by num_in_clus)
        best = max(clusters, key=lambda c: c[4])
        run_num = best[2]
        pattern = rf'DOCKED: MODEL\s+\d+\nDOCKED: USER\s+Run = {run_num}\n(.*?)(?=DOCKED: MODEL|DOCKED: ENDMDL|\Z)'
        docked_match = re.search(pattern, text, re.DOTALL)
        if not docked_match:
            return False
        lines = []
        for dl in docked_match.group(0).split('\n'):
            if dl.startswith('DOCKED: '):
                lines.append(dl[8:])
            elif dl.startswith('DOCKED:'):
                lines.append(dl[7:])
        if not lines:
            return False
        out_pdbqt.write_text('\n'.join(lines))
        return True
    except Exception:
        return False


def run_adgpu(ligand_pdbqt, maps_fld, output_dlg, output_best_pdbqt, nrun=20):
    cmd = [
        ADGPU,
        "--lfile", str(ligand_pdbqt),
        "--ffile", str(maps_fld),
        "--nrun", str(nrun),
        "--gbest", "1",
        "--xmloutput", "0",
        "--dlgoutput", "1",
    ]
    try:
        start = time.time()
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        elapsed = time.time() - start

        dlg_found = list(Path(ligand_pdbqt).parent.glob("*.dlg"))
        best_found = list(Path(ligand_pdbqt).parent.glob("*-best.pdbqt"))

        if dlg_found:
            shutil.move(str(dlg_found[0]), str(output_dlg))
        if best_found:
            shutil.move(str(best_found[0]), str(output_best_pdbqt))

        return result.returncode == 0 and output_dlg.exists(), elapsed
    except subprocess.TimeoutExpired:
        return False, 0.0
    except Exception as e:
        print(f"  ADGPU error: {e}")
        return False, 0.0


def main():
    parser = argparse.ArgumentParser(description="Run AutoDock-GPU on symmetry-corrected receptors")
    parser.add_argument("--sym-receptor-dir", required=True,
                        help="directory of symmetry-corrected PDBs (symmetry_corrected)")
    parser.add_argument("--ligand-dir", required=True, help="Meeko-prepared ligand PDBQT dir")
    parser.add_argument("--output-dir", required=True, help="output/work directory")
    parser.add_argument("--system-list", required=True)
    parser.add_argument("--system-id", type=str, default=None)
    parser.add_argument("--nrun", type=int, default=20)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    sym_dir = Path(args.sym_receptor_dir)
    ligand_dir = Path(args.ligand_dir)
    work_dir = Path(args.output_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    with open(args.system_list) as f:
        system_ids = [l.strip() for l in f if l.strip()]
    if args.system_id:
        system_ids = [s for s in system_ids if s == args.system_id]

    ok, skip, fail = 0, 0, 0
    for sys_id in system_ids:
        sym_pdb = sym_dir / sys_id / f"{sys_id}_receptor_symm.pdb"
        if not sym_pdb.exists():
            print(f"  no_receptor: {sys_id}")
            fail += 1
            continue

        lig_pdbqt = next(iter((ligand_dir / sys_id).glob("*.pdbqt")), None)
        if lig_pdbqt is None:
            print(f"  no_ligand: {sys_id}")
            fail += 1
            continue

        # Box from ligand coords
        coords = []
        with open(lig_pdbqt) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        coords.append((x, y, z))
                    except Exception:
                        continue
        if not coords:
            print(f"  no_coords: {sys_id}")
            fail += 1
            continue

        xs = [c[0] for c in coords]
        ys = [c[1] for c in coords]
        zs = [c[2] for c in coords]
        cx = (min(xs) + max(xs)) / 2.0
        cy = (min(ys) + max(ys)) / 2.0
        cz = (min(zs) + max(zs)) / 2.0
        sx = max(max(xs) - min(xs) + 10.0, 15.0)
        sy = max(max(ys) - min(ys) + 10.0, 15.0)
        sz = max(max(zs) - min(zs) + 10.0, 15.0)

        # Flat output dir <sys_id>/ — analyze script handles both flat and nested
        out_dir = work_dir / sys_id
        out_dir.mkdir(parents=True, exist_ok=True)
        out_dlg = out_dir / "output.dlg"
        out_best = out_dir / "output-best.pdbqt"
        mostpop_pdbqt = out_dir / "output-mostpop.pdbqt"
        runtime_file = out_dir / "runtime.json"

        if args.resume and out_dlg.exists() and out_dlg.stat().st_size > 0:
            skip += 1
            continue

        print(f"  Prep+dock {sys_id}...", end="", flush=True)
        tmp_dir = out_dir / "_tmp"
        if tmp_dir.exists():
            shutil.rmtree(tmp_dir)
        tmp_dir.mkdir(parents=True, exist_ok=True)

        try:
            # Step 1: generate PDBQT + GPF from sym PDB (not timed — no analog in other methods)
            rec_pdbqt = tmp_dir / "rec.pdbqt"
            gpf_path = tmp_dir / "rec.gpf"
            glg_path = tmp_dir / "rec.glg"
            rec_result = subprocess.run(
                [*MEEKO_REC.split(),
                 "--read_pdb", str(sym_pdb),
                 "-o", str(tmp_dir / "rec"),
                 "--box_center", f"{cx:.3f}", f"{cy:.3f}", f"{cz:.3f}",
                 "--box_size", f"{sx:.3f}", f"{sy:.3f}", f"{sz:.3f}",
                 "--delete_bad_res", "--default_altloc", "A",
                 "-p", "-g"],
                capture_output=True, text=True, timeout=300,
            )
            if not rec_pdbqt.exists() or not gpf_path.exists():
                print(f" meeko_fail: {rec_result.stderr[:100]}")
                fail += 1
                continue

            rec_types = get_receptor_types(rec_pdbqt)
            lig_actual_types = get_ligand_types(lig_pdbqt)
            all_types = list(dict.fromkeys(rec_types + [t for t in lig_actual_types if t not in rec_types]))
            fix_gpf(gpf_path, rec_types, all_types)

            # Step 2: autogrid4 (CPU, not timed separately)
            if not run_autogrid4(gpf_path, glg_path):
                print(" autogrid_fail")
                fail += 1
                continue

            maps_fld = gpf_path.with_name(gpf_path.stem + ".maps.fld")
            if not maps_fld.exists():
                print(" no_maps")
                fail += 1
                continue

            # Step 3: AutoDock-GPU (timed — this is the comparable docking time)
            success, elapsed = run_adgpu(
                lig_pdbqt, maps_fld, out_dlg, out_best, nrun=args.nrun,
            )
            if success:
                # Also extract the most-populated cluster's best pose
                _extract_mostpop_pose(out_dlg, mostpop_pdbqt)

                with open(runtime_file, "w") as f:
                    json.dump({"method": "autodock-gpu", "nrun": args.nrun,
                               "runtime_seconds": elapsed}, f)
                print(f" OK ({elapsed:.1f}s)")
                ok += 1
            else:
                print(" adgpu_fail")
                fail += 1
        except Exception as e:
            print(f" error: {e}")
            fail += 1
        finally:
            if tmp_dir.exists():
                shutil.rmtree(tmp_dir)

    print(f"AutoDock-GPU: OK={ok}, Skip={skip}, Fail={fail}")


if __name__ == "__main__":
    main()
