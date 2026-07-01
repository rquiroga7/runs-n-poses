#!/usr/bin/env python3
"""Count protein metal atoms within 6A of any ligand atom for systems list.

Writes a TSV with `system_id\tmetal_count` to the output file.
Relies on `obabel-25-07` being available to convert CIF/SDF -> PDB for parsing.
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path
import math

METALS = {
    'ZN','MG','CA','FE','MN','NA','K','CU','CO','NI','CD','SR','LI','AL','PB','HG',
    'TI','V','CR','MO','SE','CS','BA','AG','AU'
}


def parse_pdb_atoms(pdb_path):
    """Return list of (element, (x,y,z)) for ATOM/HETATM lines."""
    atoms = []
    try:
        with open(pdb_path, 'r') as fh:
            for line in fh:
                if not line.startswith(('ATOM', 'HETATM')):
                    continue
                if len(line) < 54:
                    continue
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                except Exception:
                    continue
                # element in columns 76-78 if present
                elem = ''
                if len(line) >= 78:
                    elem = line[76:78].strip()
                if not elem:
                    # fallback: derive from atom name (cols 12-16)
                    name = line[12:16].strip()
                    # element is letters from name
                    letters = ''.join([c for c in name if c.isalpha()])
                    elem = letters[:2].upper()
                atoms.append((elem.upper(), (x, y, z)))
    except Exception:
        pass
    return atoms


def convert_with_obabel(in_path: Path, out_path: Path, out_format: str = 'pdb'):
    cmd = ["obabel-25-07", str(in_path), "-O", str(out_path)]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        return r.returncode == 0 and out_path.exists()
    except Exception:
        return False


def distance(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)


def count_metals_for_system(system_id: str, datasets_root: Path) -> int:
    gt = datasets_root / 'ground_truth' / system_id
    cif = gt / 'receptor.cif'
    lig_dir = gt / 'ligand_files'
    if not cif.exists() or not lig_dir.exists():
        return 0

    tmp_dir = Path('/tmp') / f'rnp_metal_check_{system_id}'
    tmp_dir.mkdir(parents=True, exist_ok=True)
    pdb_rec = tmp_dir / 'rec.pdb'
    ok = convert_with_obabel(cif, pdb_rec)
    if not ok:
        return 0
    rec_atoms = parse_pdb_atoms(pdb_rec)
    # select metal atoms (keep coordinates)
    metal_atoms = [(e, coord) for e, coord in rec_atoms if e in METALS]
    if not metal_atoms:
        return 0

    # gather ligand atom coords across all ligand files
    ligand_coords = []
    for sdf in sorted(lig_dir.glob('*.sdf')):
        tmp_lig = tmp_dir / (sdf.stem + '.pdb')
        ok = convert_with_obabel(sdf, tmp_lig)
        if not ok:
            continue
        lig_atoms = parse_pdb_atoms(tmp_lig)
        ligand_coords.extend([coord for _, coord in lig_atoms])

    if not ligand_coords:
        return 0

    # count unique metal atoms within 6.0 A of any ligand atom
    count = 0
    for _, mcoord in metal_atoms:
        near = any(distance(mcoord, lcoord) <= 6.0 for lcoord in ligand_coords)
        if near:
            count += 1
    return count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--system-list', type=str, required=True)
    parser.add_argument('--datasets-root', type=str, default='/home/rquiroga/Datasets/runs-n-poses-datasets')
    parser.add_argument('--output', type=str, default='scripts/single_ligand_systems_metals.txt')
    args = parser.parse_args()

    systems = [l.strip() for l in open(args.system_list) if l.strip()]
    out_lines = []
    datasets_root = Path(args.datasets_root)
    for sid in systems:
        n = count_metals_for_system(sid, datasets_root)
        out_lines.append(f"{sid}\t{n}\n")

    outp = Path(args.output)
    outp.parent.mkdir(parents=True, exist_ok=True)
    outp.write_text(''.join(out_lines))
    print(f"Wrote {len(out_lines)} lines to {outp}")


if __name__ == '__main__':
    main()
