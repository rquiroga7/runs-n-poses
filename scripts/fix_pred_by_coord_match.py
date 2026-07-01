#!/usr/bin/env python3
"""Fix predicted ligand by mapping predicted coordinates onto GT topology.

Usage: python scripts/fix_pred_by_coord_match.py <pred_sdf|pred_pdbt> <gt_sdf> <out_sdf>

This script prefers a topology-based atom mapping from the docked SDF and
falls back to parsing coordinates from a PDBT/PDB-like file if needed.
The repaired output keeps the GT topology while transferring the docked pose
coordinates onto the matched atoms.
"""
import sys
import math
from rdkit import Chem


def is_float(s):
    try:
        float(s)
        return True
    except Exception:
        return False


def load_first_mol(path):
    suppl = Chem.SDMolSupplier(path, sanitize=False, removeHs=False)
    for m in suppl:
        if m is not None:
            return m
    return None


def strip_hydrogens_no_sanitize(mol):
    rw = Chem.RWMol(mol)
    atom_indices = [a.GetIdx() for a in rw.GetAtoms() if a.GetSymbol() == 'H']
    for idx in sorted(atom_indices, reverse=True):
        rw.RemoveAtom(idx)
    return rw.GetMol()


def parse_pdbt_coords(path):
    coords = []  # list of (element, x,y,z)
    with open(path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            # Typical layout: ATOM idx name res ... x y z ...
            parts = line.split()
            if len(parts) < 5:
                continue
            atom_name = parts[2]
            # element guess: first char of atom_name (safe fallback)
            element = atom_name[0]
            # PDB/PDBQT-style coordinates live in fixed columns, not the first
            # numeric triple in the split tokens. Using fixed widths avoids
            # accidentally reading the residue number as x.
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except Exception:
                # Fallback for slightly nonstandard spacing.
                coord_tokens = [token for token in parts if is_float(token)]
                if len(coord_tokens) < 3:
                    continue
                x, y, z = map(float, coord_tokens[:3])
            coords.append((element, x, y, z))
    return coords


def dist(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)


def map_pred_to_gt_by_graph(pred_mol, gt_mol):
    pred_heavy_indices = [a.GetIdx() for a in pred_mol.GetAtoms() if a.GetSymbol() != 'H']
    gt_heavy_indices = [a.GetIdx() for a in gt_mol.GetAtoms() if a.GetSymbol() != 'H']
    pred_heavy = strip_hydrogens_no_sanitize(pred_mol)
    gt_heavy = strip_hydrogens_no_sanitize(gt_mol)

    match = pred_heavy.GetSubstructMatch(gt_heavy)
    if match:
        return {
            gt_heavy_indices[gt_idx]: pred_heavy_indices[pred_idx]
            for gt_idx, pred_idx in enumerate(match)
        }

    reverse = gt_heavy.GetSubstructMatch(pred_heavy)
    if reverse:
        return {
            gt_heavy_indices[gt_idx]: pred_heavy_indices[pred_idx]
            for pred_idx, gt_idx in enumerate(reverse)
        }

    raise RuntimeError('Could not match heavy-atom scaffolds between pred and gt')


def map_pred_to_gt_from_coords(pred_coords, gt_mol):
    # Build lists of heavy atoms (element, idx) for pred and gt
    pred_heavy = [(i, e, x, y, z) for i, (e, x, y, z) in enumerate(pred_coords) if e.upper() != 'H']
    gt_heavy = [(a.GetIdx(), a.GetSymbol()) for a in gt_mol.GetAtoms() if a.GetSymbol() != 'H']

    mapping = {}  # gt_idx -> (x,y,z)
    used_pred = set()

    # For each GT heavy atom, find nearest pred heavy atom of same element
    for gt_idx, gt_el in gt_heavy:
        best = None
        best_d = 1e9
        for pi, pel, px, py, pz in pred_heavy:
            if pi in used_pred:
                continue
            if pel.upper() != gt_el.upper():
                continue
            d = 0.0
            try:
                conf = gt_mol.GetConformer()
                gp = conf.GetAtomPosition(gt_idx)
                d = dist((px,py,pz), (gp.x, gp.y, gp.z))
            except Exception:
                d = 0.0
            if d < best_d:
                best = (px,py,pz,pi)
                best_d = d
        if best is not None:
            px,py,pz,pi = best
            mapping[gt_idx] = (px,py,pz)
            used_pred.add(pi)

    return mapping


def assign_coords_and_write(gt_sdf, mapping, out_sdf):
    suppl = Chem.SDMolSupplier(gt_sdf, removeHs=False)
    mols = [m for m in suppl if m is not None]
    if not mols:
        raise RuntimeError('Could not read GT SDF')
    m = mols[0]

    # Ensure conformer exists
    if m.GetNumConformers() == 0:
        conf = Chem.Conformer(m.GetNumAtoms())
        for i in range(m.GetNumAtoms()):
            conf.SetAtomPosition(i, Chem.rdGeometry.Point3D(0.0,0.0,0.0))
        m.AddConformer(conf, assignId=True)

    conf = m.GetConformer()
    for gt_idx, pos in mapping.items():
        x,y,z = pos
        try:
            conf.SetAtomPosition(gt_idx, Chem.rdGeometry.Point3D(x,y,z))
        except Exception:
            pass

    w = Chem.SDWriter(out_sdf)
    w.write(m)
    w.close()


def assign_pred_mol_coords_and_write(gt_sdf, pred_mol, mapping, out_sdf):
    suppl = Chem.SDMolSupplier(gt_sdf, removeHs=False)
    mols = [m for m in suppl if m is not None]
    if not mols:
        raise RuntimeError('Could not read GT SDF')
    m = mols[0]

    if m.GetNumConformers() == 0:
        conf = Chem.Conformer(m.GetNumAtoms())
        for i in range(m.GetNumAtoms()):
            conf.SetAtomPosition(i, Chem.rdGeometry.Point3D(0.0, 0.0, 0.0))
        m.AddConformer(conf, assignId=True)

    pred_conf = pred_mol.GetConformer() if pred_mol.GetNumConformers() > 0 else None
    conf = m.GetConformer()
    for gt_idx, pred_idx in mapping.items():
        try:
            if pred_conf is not None:
                pos = pred_conf.GetAtomPosition(pred_idx)
                conf.SetAtomPosition(gt_idx, pos)
        except Exception:
            pass

    w = Chem.SDWriter(out_sdf)
    w.write(m)
    w.close()


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(2)
    pred_path, gt_path, out_path = sys.argv[1:4]
    gt_mol = load_first_mol(gt_path)
    if gt_mol is None:
        print('Failed to load GT SDF')
        sys.exit(2)

    if pred_path.lower().endswith('.sdf'):
        pred_mol = load_first_mol(pred_path)
        if pred_mol is not None:
            try:
                mapping = map_pred_to_gt_by_graph(pred_mol, gt_mol)
                assign_pred_mol_coords_and_write(gt_path, pred_mol, mapping, out_path)
                print('Wrote fixed SDF to', out_path)
                return
            except Exception:
                pass

    pred_coords = parse_pdbt_coords(pred_path)
    if not pred_coords:
        print('No coordinates parsed from predicted file')
        sys.exit(1)
    mapping = map_pred_to_gt_from_coords(pred_coords, gt_mol)
    if not mapping:
        print('Failed to map predicted coords to GT heavy atoms')
        sys.exit(3)
    assign_coords_and_write(gt_path, mapping, out_path)
    print('Wrote fixed SDF to', out_path)


if __name__ == '__main__':
    main()
