#!/usr/bin/env python3
"""Fix predicted ligand SDF by copying valence/charge/H layout from ground-truth using RDKit.

Usage: python3 scripts/fix_pred_copy_gt_valence.py <pred_sdf> <gt_sdf> <out_sdf>
"""
import sys
from rdkit import Chem


def load_first_mol(path):
    for sanitize in (False, True):
        suppl = Chem.SDMolSupplier(path, sanitize=sanitize, removeHs=False)
        for m in suppl:
            if m is not None:
                return m
    for sanitize in (False, True):
        suppl = Chem.SDMolSupplier(path, sanitize=sanitize, removeHs=True)
        for m in suppl:
            if m is not None:
                return m
    return None


def heavy_atom_index_map(mol):
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() != 'H']


def build_heavy_to_hydrogens(mol):
    d = {}
    for a in mol.GetAtoms():
        if a.GetSymbol() == 'H':
            for nb in a.GetNeighbors():
                if nb.GetSymbol() != 'H':
                    d.setdefault(nb.GetIdx(), []).append(a.GetIdx())
    return d


def copy_coords_from_pred_to_gt_topology(pred_mol, gt_mol):
    # Get heavy-atom lists and mapping indices
    pred_heavy = heavy_atom_index_map(pred_mol)
    gt_heavy = heavy_atom_index_map(gt_mol)

    pred_nohs = Chem.RemoveHs(pred_mol)
    gt_nohs = Chem.RemoveHs(gt_mol)

    match = pred_nohs.GetSubstructMatch(gt_nohs)
    if not match:
        # try reverse matching
        match = gt_nohs.GetSubstructMatch(pred_nohs)
        if not match:
            raise RuntimeError("Could not match heavy-atom scaffolds between pred and gt")
        # reverse mapping: indices in pred_nohs correspond to gt_nohs atoms
        # construct map gt_idx -> pred_idx
        inv = {pred_idx: gt_idx for gt_idx, pred_idx in enumerate(match)}
        mapping = {}
        for i, gt_orig in enumerate(gt_heavy):
            if i in inv:
                mapping[gt_orig] = pred_heavy[inv[i]]
    else:
        # match is tuple mapping gt_nohs atom i -> pred_nohs atom index
        mapping = {}
        for i, gt_orig in enumerate(gt_heavy):
            pred_nohs_idx = match[i]
            mapping[gt_orig] = pred_heavy[pred_nohs_idx]

    # Build new molecule from gt topology (copy) and assign coordinates from pred
    new = Chem.Mol(gt_mol)
    # Ensure we have a conformer to assign
    pred_conf = None
    if pred_mol.GetNumConformers() > 0:
        pred_conf = pred_mol.GetConformer()
    gt_conf = gt_mol.GetConformer() if gt_mol.GetNumConformers() > 0 else None

    new.RemoveAllConformers()
    conf = Chem.Conformer(new.GetNumAtoms())

    # copy heavy atom coords
    for gt_idx, pred_idx in mapping.items():
        try:
            if pred_conf is not None:
                pos = pred_conf.GetAtomPosition(pred_idx)
            elif gt_conf is not None:
                pos = gt_conf.GetAtomPosition(gt_idx)
            else:
                pos = Chem.rdGeometry.Point3D(0.0, 0.0, 0.0)
            conf.SetAtomPosition(gt_idx, pos)
        except Exception:
            pass

    # For hydrogens, prefer positions from GT if available, otherwise leave at 0
    gt_h_map = build_heavy_to_hydrogens(gt_mol)
    pred_h_map = build_heavy_to_hydrogens(pred_mol)
    for heavy_idx, h_list in gt_h_map.items():
        for i, h_idx in enumerate(h_list):
            # try to copy from GT conf
            if gt_conf is not None:
                try:
                    pos = gt_conf.GetAtomPosition(h_idx)
                    conf.SetAtomPosition(h_idx, pos)
                    continue
                except Exception:
                    pass
            # fallback: if pred has H attached to corresponding pred heavy
            pred_heavy_idx = mapping.get(heavy_idx)
            if pred_heavy_idx is not None and pred_h_map.get(pred_heavy_idx):
                try:
                    pred_h_idx = pred_h_map[pred_heavy_idx][0]
                    if pred_conf is not None:
                        pos = pred_conf.GetAtomPosition(pred_h_idx)
                        conf.SetAtomPosition(h_idx, pos)
                        continue
                except Exception:
                    pass
            # last resort: copy heavy atom pos
            try:
                heavy_pos = conf.GetAtomPosition(heavy_idx)
                conf.SetAtomPosition(h_idx, heavy_pos)
            except Exception:
                pass

    new.AddConformer(conf, assignId=True)
    return new


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(2)
    pred_path, gt_path, out_path = sys.argv[1:4]
    pred = load_first_mol(pred_path)
    gt = load_first_mol(gt_path)
    if pred is None or gt is None:
        print('Failed to load pred or gt')
        sys.exit(1)
    try:
        fixed = copy_coords_from_pred_to_gt_topology(pred, gt)
    except Exception as e:
        print('Failed to reconcile:', e)
        sys.exit(2)

    # write out
    w = Chem.SDWriter(out_path)
    w.write(fixed)
    w.close()


if __name__ == '__main__':
    main()
