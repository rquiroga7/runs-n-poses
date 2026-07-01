#!/usr/bin/env python
"""
Run AutoDock-GPU (Scripps) on symmetry-corrected receptors.

Pipeline per system:
  1. Meeko mk_prepare_ligand → ligand PDBQT
  2. Meeko mk_prepare_receptor → receptor PDBQT + GPF
  3. Fix GPF: filter receptor_types to actual PDBQT types, drop unused ligand types
  4. autogrid4 → grid maps (.maps.fld)
  5. AutoDock-GPU → docking results (.dlg)

Usage:
    python run_autodock_gpu_symmetry.py \
        --sym-dir runs-n-poses-datasets/symmetry_corrected \
        --work-dir runs-n-poses-datasets/autodock_gpu_symmetry \
        --system-list scripts/systems_for_symmetry_docking.txt
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

ADGPU = "AutoDock-GPU"
AUTOGRID4 = "/home/rquiroga/Downloads/autodocksuite-4.2.6-x86_64Linux2/x86_64Linux2/autogrid4"
MEEKOREC = "python3 -m meeko.cli.mk_prepare_receptor"
MEEKOLIG = "python3 -m meeko.cli.mk_prepare_ligand"
PARAM_FILE = "boron-silicon-atom_par.dat"


def get_receptor_types(pdbqt_path):
    """Extract unique AD4 atom types from a receptor PDBQT."""
    types = set()
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                parts = line.split()
                if len(parts) >= 12:
                    types.add(parts[-1])
    return sorted(types)


def get_ligand_types(pdbqt_path):
    """Extract unique AD4 atom types from a ligand PDBQT."""
    types = set()
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                parts = line.split()
                if len(parts) >= 12:
                    types.add(parts[-1])
    return sorted(types)


def fix_gpf(gpf_path, rec_types, lig_types):
    """Fix GPF so autogrid4 accepts it.

    autogrid4 v4.2.6 requires:
      - receptor_types to only list types present in the receptor PDBQT
      - ligand_types + map count to match exactly
    """
    with open(gpf_path) as f:
        lines = f.readlines()

    new_lines = []
    seen_types = False
    for line in lines:
        if line.startswith("receptor_types"):
            new_lines.append("receptor_types " + " ".join(rec_types) + "\n")
            seen_types = True
        else:
            new_lines.append(line)

    assert seen_types, "No receptor_types line found in GPF"

    gpf = "".join(new_lines)

    # Filter map lines to only types present in ligand
    # Also ensure map count matches ligand_types
    gpf = re.sub(
        r"^ligand_types .*",
        "ligand_types " + " ".join(lig_types),
        gpf,
        flags=re.MULTILINE,
    )

    # Keep only map/elecmap/dsolvmap lines for types we actually have
    valid_maps = set(lig_types) | {"e", "d"}
    def keep_map(m):
        type_ = m.group(1)
        if type_ in valid_maps:
            return m.group(0)
        return ""

    gpf = re.sub(r"^map .*?\.(\w+)\.map", lambda m: m.group(0) if m.group(1) in lig_types else "", gpf, flags=re.MULTILINE)
    gpf = re.sub(r"\n{2,}", "\n", gpf)
    gpf = gpf.strip() + "\n"

    with open(gpf_path, "w") as f:
        f.write(gpf)


def run_autogrid4(gpf_path, glg_path):
    """Run autogrid4 to generate grid maps."""
    result = subprocess.run(
        [AUTOGRID4, "-p", str(gpf_path), "-l", str(glg_path)],
        capture_output=True, text=True, timeout=300,
    )
    return result.returncode == 0


def run_adgpu(ligand_pdbqt, maps_fld, output_dlg, output_best_pdbqt, nrun=20):
    """Run AutoDock-GPU docking."""
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

        # AutoDock-GPU writes output files to same dir as ligand
        dlg_found = list(Path(ligand_pdbqt).parent.glob("*.dlg"))
        best_found = list(Path(ligand_pdbqt).parent.glob("*-best.pdbqt"))

        if dlg_found:
            shutil.move(str(dlg_found[0]), str(output_dlg))
        if best_found:
            shutil.move(str(best_found[0]), str(output_best_pdbqt))

        return result.returncode == 0 and output_dlg.exists(), elapsed
    except subprocess.TimeoutExpired:
        return False, 0
    except Exception as e:
        print(f"  ADGPU error: {e}")
        return False, 0


def main():
    parser = argparse.ArgumentParser(description="Run AutoDock-GPU on symmetry-corrected receptors")
    parser.add_argument("--sym-dir", required=True, help="symmetry_corrected directory")
    parser.add_argument("--work-dir", required=True, help="output/work directory")
    parser.add_argument("--system-list", required=True)
    parser.add_argument("--system-id", type=str, default=None)
    parser.add_argument("--nrun", type=int, default=20)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    sym_dir = Path(args.sym_dir)
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    with open(args.system_list) as f:
        system_ids = [l.strip() for l in f if l.strip()]
    if args.system_id:
        system_ids = [s for s in system_ids if s == args.system_id]

    ok, fail, skip = 0, 0, 0
    for sys_id in system_ids:
        sys_dir = sym_dir / sys_id
        if not sys_dir.exists():
            print(f"  No system dir: {sys_id}")
            fail += 1
            continue

        ligand_pdb = list(sys_dir.glob("*_ligand.pdb"))
        receptor_pdb = list(sys_dir.glob("*_receptor_symm.pdb"))
        if not ligand_pdb or not receptor_pdb:
            print(f"  Missing input files: {sys_id}")
            fail += 1
            continue

        out_dir = work_dir / sys_id
        out_dir.mkdir(parents=True, exist_ok=True)
        out_dlg = out_dir / "output.dlg"
        out_best = out_dir / "output-best.pdbqt"
        runtime_file = out_dir / "runtime.json"

        if args.resume and out_dlg.exists() and out_dlg.stat().st_size > 0:
            skip += 1
            continue

        tmp_dir = out_dir / "_tmp"
        tmp_dir.mkdir(parents=True, exist_ok=True)

        try:
            # Step 1: Meeko ligand prep
            lig_pdbqt = tmp_dir / "ligand.pdbqt"
            lig_result = subprocess.run(
                [*MEEKOLIG.split(), "-i", str(ligand_pdb[0]), "-o", str(lig_pdbqt)],
                capture_output=True, text=True, timeout=120,
            )
            if not lig_pdbqt.exists():
                print(f"  Ligand prep failed: {sys_id}")
                print(f"    {lig_result.stderr[:200]}")
                fail += 1
                continue

            # Step 2: Meeko receptor prep (get box from ligand coordinates)
            coords = []
            with open(ligand_pdb[0]) as f:
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
                print(f"  No coordinates: {sys_id}")
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

            gpf_path = tmp_dir / "receptor.gpf"
            glg_path = tmp_dir / "receptor.glg"
            rec_pdbqt = tmp_dir / "receptor.pdbqt"

            rec_result = subprocess.run(
                [*MEEKOREC.split(),
                 "-i", str(receptor_pdb[0]),
                 "-o", str(tmp_dir / "rec"),
                 f"--box_center={cx:.3f},{cy:.3f},{cz:.3f}",
                 f"--box_size={sx:.3f},{sy:.3f},{sz:.3f}",
                 "--delete_bad_res", "--default_altloc", "A",
                 "-p", "-g"],
                capture_output=True, text=True, timeout=300,
            )
            if not rec_pdbqt.exists() or not gpf_path.exists():
                print(f"  Receptor prep failed: {sys_id}")
                print(f"    {rec_result.stderr[:200]}")
                fail += 1
                continue

            # Step 3: Fix GPF
            rec_types = get_receptor_types(rec_pdbqt)
            lig_actual_types = get_ligand_types(lig_pdbqt)
            # Merge rec + lig types for ligand_types; remove any that exceed AD4 limit
            all_types = list(dict.fromkeys(rec_types + [t for t in lig_actual_types if t not in rec_types]))
            # autogrid4 seems to have a 14-type limit for ligand_types
            # Keep only types that are actually present in either receptor or ligand
            fix_gpf(gpf_path, rec_types, all_types)

            # Step 4: autogrid4
            if not run_autogrid4(gpf_path, glg_path):
                print(f"  autogrid4 failed: {sys_id}")
                fail += 1
                continue

            maps_fld = gpf_path.with_name(gpf_path.stem + ".maps.fld")
            if not maps_fld.exists():
                print(f"  No maps.fld: {sys_id}")
                fail += 1
                continue

            # Step 5: AutoDock-GPU
            success, elapsed = run_adgpu(
                lig_pdbqt, maps_fld, out_dlg, out_best, nrun=args.nrun,
            )

            if success:
                with open(runtime_file, "w") as f:
                    json.dump({"method": "autodock-gpu", "nrun": args.nrun,
                               "runtime_seconds": elapsed}, f)
                ok += 1
            else:
                print(f"  ADGPU failed: {sys_id}")
                fail += 1

        except Exception as e:
            print(f"  Error {sys_id}: {e}")
            fail += 1
        finally:
            if tmp_dir.exists():
                shutil.rmtree(tmp_dir)

    print(f"AutoDock-GPU: OK={ok}, Skip={skip}, Fail={fail}")


if __name__ == "__main__":
    main()
