#!/usr/bin/env python
"""
Step 2b: Meeko-based ligand preparation for runs-n-poses.

Reads the proper ligand SDF from ``ground_truth/<sys_id>/ligand_files/<chain>.sdf``
and converts to both PDBT (for Vinardo) and PDBQT (for Vina) using Meeko's
``mk_prepare_pdbt_ligand`` and ``mk_prepare_ligand`` CLIs. The PDBT writer
is part of the new Meeko release at ``/home/rquiroga/github/Meeko``.

Meeko provides more accurate atom typing and macrocycle/torsion handling
than the obabel-based pipeline used in ``prepare_symmetry_ligands.py``.

The Meeko calls are expensive (a few seconds per system), so we run
them in parallel via ``concurrent.futures.ProcessPoolExecutor``.

Usage:
    python prepare_ligand_meeko.py \
        --ground-truth-dir runs-n-poses-datasets/ground_truth \
        --output-pdbt runs-n-poses-datasets/meeko_ligands_pdbt \
        --output-pdbqt runs-n-poses-datasets/meeko_ligands_pdbqt \
        --system-list scripts/single_ligand_systems_symmetry.csv
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from rdkit import Chem

MEEKO_PDBT = "/home/rquiroga/github/runs-n-poses/.venv/bin/python -m meeko.cli.mk_prepare_pdbt_ligand"
MEEKO_PDBQT = "/home/rquiroga/github/runs-n-poses/.venv/bin/python -m meeko.cli.mk_prepare_ligand"


def _add_hydrogens_sdf(sdf_path: str) -> str:
    """Ensure the SDF has explicit Hs (Meeko requires them)."""
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        return sdf_path
    if any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
        return sdf_path
    mol = Chem.AddHs(mol, addCoords=True)
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
    tmp.close()
    writer = Chem.SDWriter(tmp.name)
    writer.write(mol)
    writer.close()
    return tmp.name


def _run(meeko_cmd: str, input_sdf: str, output_path: str) -> bool:
    cmd = meeko_cmd.split() + ["-i", input_sdf, "-o", output_path]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        return r.returncode == 0 and os.path.exists(output_path) and os.path.getsize(output_path) > 0
    except Exception:
        return False


def _convert_to_pdbt(sdf_path: str, output_path: str) -> bool:
    h_sdf = _add_hydrogens_sdf(sdf_path)
    try:
        return _run(MEEKO_PDBT, h_sdf, output_path)
    finally:
        if h_sdf != sdf_path and os.path.exists(h_sdf):
            os.remove(h_sdf)


def _convert_to_pdbqt(sdf_path: str, output_path: str) -> bool:
    h_sdf = _add_hydrogens_sdf(sdf_path)
    try:
        return _run(MEEKO_PDBQT, h_sdf, output_path)
    finally:
        if h_sdf != sdf_path and os.path.exists(h_sdf):
            os.remove(h_sdf)


def _worker(args):
    """Module-level worker (picklable) — runs both PDBT and PDBQT prep."""
    sys_id, chain, gt_dir, out_pdbt, out_pdbqt = args
    gt_dir = Path(gt_dir)
    out_pdbt = Path(out_pdbt)
    out_pdbqt = Path(out_pdbqt)

    sdf_path = gt_dir / sys_id / "ligand_files" / f"{chain}.sdf"
    if not sdf_path.exists():
        return (sys_id, chain), False, False, "no_sdf"

    pdbt_dir = out_pdbt / sys_id
    pdbt_dir.mkdir(parents=True, exist_ok=True)
    pdbt_path = pdbt_dir / f"{chain}.pdbt"
    if not (pdbt_path.exists() and pdbt_path.stat().st_size > 0):
        ok1 = _convert_to_pdbt(str(sdf_path), str(pdbt_path))
    else:
        ok1 = True

    pdbqt_dir = out_pdbqt / sys_id
    pdbqt_dir.mkdir(parents=True, exist_ok=True)
    pdbqt_path = pdbqt_dir / f"{chain}.pdbqt"
    if not (pdbqt_path.exists() and pdbqt_path.stat().st_size > 0):
        ok2 = _convert_to_pdbqt(str(sdf_path), str(pdbqt_path))
    else:
        ok2 = True

    return (sys_id, chain), ok1, ok2, "ok"


def main():
    parser = argparse.ArgumentParser(description="Meeko-based ligand preparation")
    parser.add_argument("--ground-truth-dir", required=True)
    parser.add_argument("--output-pdbt", required=True)
    parser.add_argument("--output-pdbqt", required=True)
    parser.add_argument("--system-info", required=True,
                        help="CSV with system_id,proper_ligand_chain (matches prepare_symmetry_ligands.py)")
    parser.add_argument("--system-id", type=str, default=None)
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4)))
    args = parser.parse_args()

    gt_dir = Path(args.ground_truth_dir)
    out_pdbt = Path(args.output_pdbt)
    out_pdbqt = Path(args.output_pdbqt)
    out_pdbt.mkdir(parents=True, exist_ok=True)
    out_pdbqt.mkdir(parents=True, exist_ok=True)

    with open(args.system_info) as f:
        reader = csv.DictReader(f)
        systems = [(r["system_id"], r["proper_ligand_chain"]) for r in reader]
    if args.system_id:
        systems = [s for s in systems if s[0] == args.system_id]

    try:
        from tqdm import tqdm
        progress = tqdm(total=len(systems), desc="Meeko ligand")
    except ImportError:
        progress = None

    ok_pdbt, ok_pdbqt, fail = 0, 0, 0
    args_iter = [(s[0], s[1], str(gt_dir), str(out_pdbt), str(out_pdbqt)) for s in systems]
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futures = {ex.submit(_worker, a): a[:2] for a in args_iter}
        for fut in as_completed(futures):
            try:
                item, ok1, ok2, status = fut.result()
            except Exception as e:
                if progress:
                    progress.update(1)
                fail += 1
                continue
            if status != "ok":
                fail += 1
            else:
                if ok1:
                    ok_pdbt += 1
                else:
                    print(f"  FAIL pdbt: {item[0]} {item[1]}")
                    fail += 1
                if ok2:
                    ok_pdbqt += 1
                else:
                    print(f"  FAIL pdbqt: {item[0]} {item[1]}")
                    fail += 1
            if progress:
                progress.update(1)
    if progress:
        progress.close()

    print(f"Ligands (Meeko): PDBT={ok_pdbt}, PDBQT={ok_pdbqt}, Fail={fail}")


if __name__ == "__main__":
    main()
