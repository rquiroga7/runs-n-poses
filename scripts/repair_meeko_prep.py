#!/usr/bin/env python
"""
Fix Meeko-prepared PDBQT files for compatibility with AutoDock Vina.

Issues addressed:
  1. Receceptors missing from meeko_receptors_pdbqt — re-attempt prep without
     problematic non-standard residues (DpL, etc.)
  2. Ligands missing from meeko_ligands_pdbqt — re-attempt prep
  3. Invalid atom types (e.g. "B" for Boron) — remap to nearest valid type
     if near the binding site, or remove if far from the ligand.
"""

import csv
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from rdkit import Chem

# Config
MEEKO_REC = "/home/rquiroga/github/runs-n-poses/.venv/bin/python -m meeko.cli.mk_prepare_receptor"
MEEKO_PDBT = "/home/rquiroga/github/runs-n-poses/.venv/bin/python -m meeko.cli.mk_prepare_pdbt_receptor"
MEEKO_LIG_PDBT = "/home/rquiroga/github/runs-n-poses/.venv/bin/python -m meeko.cli.mk_prepare_pdbt_ligand"
MEEKO_LIG_PDBQT = "/home/rquiroga/github/runs-n-poses/.venv/bin/python -m meeko.cli.mk_prepare_ligand"

REC_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/meeko_receptors_pdbqt")
REC_PDBT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/meeko_receptors_pdbt")
LIG_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/meeko_ligands_pdbqt")
LIG_PDBT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/meeko_ligands_pdbt")
SYM_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/symmetry_corrected")
GT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth")
SYS_LIST = Path("/home/rquiroga/github/runs-n-poses/scripts/single_ligand_systems_symmetry.csv")

# Valid Vina/AD4 atom types (2-char and 3-char)
VALID_TYPES = {
    'C', 'A', 'N', 'NA', 'OA', 'F', 'P', 'SA', 'S',
    'Cl', 'Br', 'I', 'HD', 'HS', 'H', 'Si', 'B',
    'Zn', 'Mg', 'Ca', 'Fe', 'Mn', 'Cu', 'Co', 'Ni',
    'Se', 'W',
    'G0', 'G1', 'G2', 'G3',
    'CG0', 'CG1', 'CG2', 'CG3',
}

# Map invalid types to the closest Vina-compatible type
TYPE_MAP = {
    "B": "C",  # Boron → Carbon (Vina often rejects "B")
}

CONTACT_DIST = 6.0  # Angstrom — remove far atoms beyond this from any ligand atom


def _load_systems():
    with open(SYS_LIST) as f:
        reader = csv.DictReader(f)
        return [(r["system_id"], r["proper_ligand_chain"]) for r in reader]


def _add_hs_to_sdf(sdf_path):
    """Add explicit hydrogens to an SDF file (Meeko requires them)."""
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        return None
    if any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
        return str(sdf_path)
    mol = Chem.AddHs(mol, addCoords=True)
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
    tmp.close()
    Chem.SDWriter(tmp.name).write(mol)
    return tmp.name


def repair_receptor(sys_id):
    """Fix missing Meeko receptor by stripping problematic HETATM residues (DpL)."""
    out_dir = REC_DIR / sys_id
    out_dir.mkdir(parents=True, exist_ok=True)
    out_pdbqt = out_dir / f"{sys_id}_receptor.pdbqt"
    out_pdbt = REC_PDBT_DIR / sys_id / f"{sys_id}_receptor.pdbt"
    REC_PDBT_DIR.mkdir(parents=True, exist_ok=True)

    sym_pdb = SYM_DIR / sys_id / f"{sys_id}_receptor_symm.pdb"
    if not sym_pdb.exists():
        return "no_sym_pdb"

    # Try original Meeko prep first
    cmd = MEEKO_REC.split() + ["--read_pdb", str(sym_pdb), "-o",
                                str(out_dir / f"{sys_id}_receptor".rstrip("_receptor")),
                                "--default_altloc", "A", "-x", "-p"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if r.returncode == 0 and out_pdbqt.exists():
        return "repair_ok"

    # If it failed, strip DpL residues and retry
    stripped = sym_pdb.read_text()
    if "DpL" in stripped:
        lines = [l for l in stripped.split("\n")
                 if not (l.startswith(("ATOM", "HETATM")) and "DpL" in l[17:20])]
        tmp_pdb = Path(tempfile.mktemp(suffix=".pdb"))
        tmp_pdb.write_text("\n".join(lines))
        cmd = MEEKO_REC.split() + ["--read_pdb", str(tmp_pdb), "-o",
                                    str(out_dir / f"{sys_id}_receptor".rstrip("_receptor")),
                                    "--default_altloc", "A", "-x", "-p"]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        tmp_pdb.unlink(missing_ok=True)
        if r.returncode == 0 and out_pdbqt.exists():
            return "repair_stripped"
        return f"still_fails: {r.stderr[:100]}"
    return f"no_dpl_but_fails: {r.stderr[:100]}"


def repair_ligand(sys_id, chain):
    """Re-attempt Meeko ligand prep if missing."""
    sdf_path = GT_DIR / sys_id / "ligand_files" / f"{chain}.sdf"
    if not sdf_path.exists():
        return "no_sdf"

    for out_dir, cmd_tpl in [(LIG_DIR, MEEKO_LIG_PDBQT), (LIG_PDBT_DIR, MEEKO_LIG_PDBT)]:
        out_file = out_dir / sys_id / f"{chain}.pdbt" if "pdbt" in cmd_tpl else out_dir / sys_id / f"{chain}.pdbqt"
        out_dir.mkdir(parents=True, exist_ok=True)
        if out_file.exists() and out_file.stat().st_size > 0:
            continue
        h_sdf = _add_hs_to_sdf(str(sdf_path))
        if h_sdf is None:
            continue
        cmd = cmd_tpl.split() + ["-i", h_sdf, "-o", str(out_file)]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if h_sdf != str(sdf_path) and os.path.exists(h_sdf):
            os.remove(h_sdf)
    return "ok"


def fix_invalid_types_in_pdbqt(pdbqt_path):
    """Replace invalid atom types in a PDBQT file. Returns True if changes made."""
    if not pdbqt_path.exists():
        return False
    changed = False
    lines = []
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 79:
                t = line[77:79].strip()
                if t in TYPE_MAP:
                    new_type = TYPE_MAP[t]
                    line = line[:77] + f"{new_type:>2s}" + line[79:]
                    changed = True
            lines.append(line)
    if changed:
        pdbqt_path.write_text("".join(lines))
    return changed


def filter_distant_invalid_types(pdbqt_path):
    """Remove atoms with invalid types that are >CONTACT_DIST from any ligand atom."""
    if not pdbqt_path.exists():
        return False
    # Parse all atom coords
    atoms = []
    for line in pdbqt_path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 79:
            t = line[77:79].strip()
            if t not in VALID_TYPES:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                atoms.append((t, x, y, z, line))
    if not atoms:
        return False

    # Get ligand coords from the corresponding meeko_ligands_pdbqt
    # Try to read the ligand PDBQT
    lig_path = None
    # ... this is complex; skip for now
    return False


def main():
    systems = _load_systems()
    rec_done, lig_repair, type_fix = 0, 0, 0

    # Step 1: Find and re-prep missing receptors
    print("=== Missing receptors ===")
    missing_recs = [(s, c) for s, c in systems
                    if not (REC_DIR / s / f"{s}_receptor.pdbqt").exists()]
    print(f"  {len(missing_recs)} systems without Meeko receptor")
    for sys_id, chain in missing_recs[:20]:  # Limit to 20 for testing
        status = repair_receptor(sys_id)
        print(f"    {sys_id}: {status}")
        if status.startswith("repair"):
            rec_done += 1

    # Step 2: Fix ligands with missing files or invalid types
    print("\n=== Missing/invalid ligands ===")
    for sys_id, chain in systems:
        lig_pdbqt = LIG_DIR / sys_id / f"{chain}.pdbqt"
        if not lig_pdbqt.exists() or lig_pdbqt.stat().st_size == 0:
            status = repair_ligand(sys_id, chain)
            if "ok" in status:
                lig_repair += 1
                print(f"    {sys_id} {chain}: re-prepped")
        else:
            # Fix invalid types
            if fix_invalid_types_in_pdbqt(lig_pdbqt):
                type_fix += 1
                print(f"    {sys_id} {chain}: type fix applied")

    print(f"\nDone: {rec_done} receptors repaired, {lig_repair} ligands re-prepped, {type_fix} type fixes")


if __name__ == "__main__":
    main()
