#!/usr/bin/env python
"""
Step 1b: Meeko-based receptor preparation.

Reads symm-corrected PDBs and produces PDBQT + PDBT via Meeko.
The Meeko Polymer class has been patched to gracefully skip DpL
residues instead of crashing (see polymer.py lines 1253-1255, 2002-2006).
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

MEEKO_PDBQT_REC = "/home/rquiroga/github/runs-n-poses/.venv/bin/python -m meeko.cli.mk_prepare_receptor"
MEEKO_PDBT = "/home/rquiroga/github/runs-n-poses/.venv/bin/python -m meeko.cli.mk_prepare_pdbt_receptor"


def _run(cmd: list, timeout: int = 300) -> bool:
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return r.returncode == 0
    except Exception:
        return False


def _worker(args):
    sys_id, input_dir, out_pdbqt, out_pdbt = args
    input_dir = Path(input_dir)
    out_pdbqt = Path(out_pdbqt)
    out_pdbt = Path(out_pdbt)

    pdb_path = input_dir / sys_id / f"{sys_id}_receptor_symm.pdb"
    if not pdb_path.exists():
        return sys_id, False, False, "no_pdb"

    pdbqt_dir = out_pdbqt / sys_id
    pdbqt_dir.mkdir(parents=True, exist_ok=True)
    pdbqt_path = pdbqt_dir / f"{sys_id}_receptor.pdbqt"
    ok1 = False
    if not (pdbqt_path.exists() and pdbqt_path.stat().st_size > 0):
        cmd = MEEKO_PDBQT_REC.split() + [
            "--read_pdb", str(pdb_path),
            "-o", str(pdbqt_path.with_suffix("")),
            "--forgive_extra_bonds",
            "--default_altloc", "A",
            "-x",
            "-p",
        ]
        _run(cmd)
        ok1 = pdbqt_path.exists() and pdbqt_path.stat().st_size > 0
    else:
        ok1 = True

    pdbt_dir = out_pdbt / sys_id
    pdbt_dir.mkdir(parents=True, exist_ok=True)
    pdbt_path = pdbt_dir / f"{sys_id}_receptor.pdbt"
    ok2 = False
    if not (pdbt_path.exists() and pdbt_path.stat().st_size > 0):
        cmd = MEEKO_PDBT.split() + [
            "-i", str(pdb_path),
            "-o", str(pdbt_path),
            "--forgive_extra_bonds",
            "--default_altloc", "A",
            "-x",
        ]
        _run(cmd)
        ok2 = pdbt_path.exists() and pdbt_path.stat().st_size > 0
    else:
        ok2 = True

    return sys_id, ok1, ok2, "ok"


def main():
    parser = argparse.ArgumentParser(description="Meeko-based receptor preparation")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-pdbqt", required=True)
    parser.add_argument("--output-pdbt", required=True)
    parser.add_argument("--system-list", required=True)
    parser.add_argument("--system-id", type=str, default=None)
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4)))
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    out_pdbqt = Path(args.output_pdbqt)
    out_pdbt = Path(args.output_pdbt)
    out_pdbqt.mkdir(parents=True, exist_ok=True)
    out_pdbt.mkdir(parents=True, exist_ok=True)

    with open(args.system_list) as f:
        reader = csv.DictReader(f)
        systems = [r["system_id"] for r in reader]
    if args.system_id:
        systems = [s for s in systems if s == args.system_id]

    try:
        from tqdm import tqdm
        progress = tqdm(total=len(systems), desc="Meeko receptor")
    except ImportError:
        progress = None

    ok_pdbqt, ok_pdbt, fail = 0, 0, 0
    args_iter = [(s, str(input_dir), str(out_pdbqt), str(out_pdbt)) for s in systems]
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futures = {ex.submit(_worker, a): a[0] for a in args_iter}
        for fut in as_completed(futures):
            try:
                sid, ok1, ok2, status = fut.result()
            except Exception:
                if progress:
                    progress.update(1)
                fail += 1
                continue
            if ok1:
                ok_pdbqt += 1
            else:
                fail += 1
            if ok2:
                ok_pdbt += 1
            if progress:
                progress.update(1)
    if progress:
        progress.close()

    print(f"Receptors (Meeko): PDBQT={ok_pdbqt}, PDBT={ok_pdbt}, Fail={fail}")


if __name__ == "__main__":
    main()
