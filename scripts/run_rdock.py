#!/usr/bin/env python
"""
Run rDock on symmetry-corrected receptors.

Pipeline per system:
  1. Convert receptor PDB -> MOL2
  2. Copy ligand SDF to RBT_ROOT/data/ligands/ref.sd
  3. Run rbcavity -W to generate .as cavity file
  4. Run rbdock with -i (ligand), -o (output), -r (rec pmr), -p (protocol prm)
  5. Score with OST

Usage:
    python run_rdock.py \
        --output-dir runs-n-poses-datasets/rdock \
        --system-list scripts/single_ligand_systems_symmetry.csv
"""

import argparse
import csv
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

RBT_ROOT = Path("/home/rquiroga/github/rDock")
RBDOCK = RBT_ROOT / "bin" / "rbdock"
RBCAVITY = RBT_ROOT / "bin" / "rbcavity"
SYM_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/symmetry_corrected")
GT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth")


def write_cavity_in(cav_in: Path, receptor_mol2: Path, ref_sdf: Path) -> None:
    cav_in.write_text(f"""RBT_PARAMETER_FILE_V1.00
RECEPTOR_FILE {receptor_mol2}
REFERENCE_LIGAND_FILE {ref_sdf}
CAVITY_RADIUS 6.0
SMALL_SPHERE 1.0
LARGE_SPHERE 4.0
GRID_STEP 0.5
CAVITY_MAPPING 1
BORDER 2.0
SECTION MAPPER
    SITE_MAPPER RbtLigandSiteMapper
END_SECTION
""")


def write_rec_prm(prm_path: Path, receptor_mol2: Path, cav_as: Path) -> None:
    prm_path.write_text(f"""RBT_PARAMETER_FILE_V1.00
RECEPTOR_FILE {receptor_mol2}
SECTION CAVITY
    CAVITY_FILE {cav_as}
END_SECTION
""")


def write_proto_prm(prm_path: Path, n_runs: int = 100) -> None:
    prm_path.write_text(f"""RBT_PARAMETER_FILE_V1.00
SECTION SCORE
    INTER RbtInterIdxSF.prm
    INTRA RbtIntraSF.prm
END_SECTION
SECTION GA
    POP_SIZE 50
    NUM_GA 3
    NUM_SAVE {n_runs}
END_SECTION
""")


def run_rdock(sys_id: str, chain: str, output_dir: Path, resume: bool) -> tuple[str, float | None]:
    out_dir = output_dir / sys_id / f"{sys_id}_{chain}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_sd = out_dir / "out.sd"
    log_file = out_dir / "log.txt"
    runtime_file = out_dir / "runtime.json"

    if resume and out_sd.exists() and out_sd.stat().st_size > 0:
        return "skip", None

    sym_pdb = SYM_DIR / sys_id / f"{sys_id}_receptor_symm.pdb"
    if not sym_pdb.exists():
        return "no_receptor", None

    lig_sdf = GT_DIR / sys_id / "ligand_files" / f"{chain}.sdf"
    if not lig_sdf.exists():
        return "no_ligand", None

    env = os.environ.copy()
    env["RBT_ROOT"] = str(RBT_ROOT)
    env["LD_LIBRARY_PATH"] = f"{RBT_ROOT}/lib:" + env.get("LD_LIBRARY_PATH", "")

    tmp_dir = out_dir / "_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    try:
        # Step 1: Convert receptor PDB -> MOL2
        rec_mol2 = tmp_dir / "receptor.mol2"
        r = subprocess.run(["obabel-25-07", str(sym_pdb), "-O", str(rec_mol2)],
                           capture_output=True, text=True, timeout=60)
        if r.returncode != 0 or not rec_mol2.exists():
            return "mol2_fail", None

        # Step 2: Prepare reference ligand for rbcavity
        ref_sd = out_dir / "ref.sd"
        shutil.copy2(lig_sdf, ref_sd)

        # Step 3: Generate cavity .as file (rbcavity -W)
        cav_in = tmp_dir / "cavity.in"
        write_cavity_in(cav_in, rec_mol2, ref_sd)
        r = subprocess.run([str(RBCAVITY), "-W", "-r", cav_in.name],
                           capture_output=True, text=True, timeout=120, env=env,
                           cwd=str(cav_in.parent))
        if r.returncode != 0:
            return f"cavity_fail: {(r.stderr or '')[:200]}", None

        # .as file is named <cavity_in_stem>.as in the cwd
        cav_as = cav_in.parent / "cavity.as"
        if not cav_as.exists():
            return "cavity_as_missing", None

        # Copy .as to output base dir; rbdock looks for <stem>.as where stem = sys_id split at first "."
        stem = sys_id.split(".")[0] if "." in sys_id else sys_id
        cav_as_base = output_dir / f"{stem}.as"
        shutil.copy2(cav_as, cav_as_base)

        # Step 4: Run docking
        rec_prm = tmp_dir / "rec.prm"
        proto_prm = tmp_dir / "proto.prm"
        write_rec_prm(rec_prm, rec_mol2, cav_as_base)
        write_proto_prm(proto_prm)

        start = time.time()
        r = subprocess.run([str(RBDOCK), "-i", str(lig_sdf), "-o", str(out_sd.with_suffix("")),
                            "-r", str(rec_prm), "-p", str(proto_prm), "-n", "100"],
                           capture_output=True, text=True, timeout=600, env=env)
        elapsed = time.time() - start

        with open(log_file, "w") as f:
            if r.stdout:
                f.write(r.stdout)
            if r.stderr:
                f.write("\nSTDERR:\n" + r.stderr)

        if r.returncode == 0 and out_sd.exists() and out_sd.stat().st_size > 0:
            with open(runtime_file, "w") as f:
                json.dump({"method": "rdock", "runtime_seconds": elapsed}, f)
            return "ok", elapsed
        return "dock_fail", None

    except subprocess.TimeoutExpired:
        return "timeout", None
    except Exception as e:
        return f"error: {e}", None
    finally:
        if tmp_dir.exists():
            shutil.rmtree(tmp_dir)


def main():
    parser = argparse.ArgumentParser(description="Run rDock on symmetry-corrected receptors")
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
        status, rt = run_rdock(sys_id, chain, output_dir, args.resume)
        if status == "ok":
            print(f"  OK {sys_id} ({rt:.1f}s)")
            ok += 1
        elif status == "skip":
            skip += 1
        else:
            print(f"  {status}: {sys_id}")
            fail += 1

    print(f"rDock: OK={ok}, Skip={skip}, Fail={fail}")


if __name__ == "__main__":
    main()
