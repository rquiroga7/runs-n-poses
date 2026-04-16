#!/usr/bin/env python3
# scripts/pb_rdkit_checks.py
import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def load_and_sanitize(path):
    m = Chem.MolFromMolFile(path, sanitize=False, removeHs=False)
    if m is None:
        block = open(path).read()
        m = Chem.MolFromMolBlock(block, sanitize=False, removeHs=False)
    if m is None:
        return None
    try:
        Chem.SanitizeMol(m)
    except Exception:
        try:
            mh = Chem.AddHs(m)
            Chem.SanitizeMol(mh)
            m = Chem.RemoveHs(mh)
        except Exception:
            return None
    return m

def mol_formula(m):
    return rdMolDescriptors.CalcMolFormula(m)

def bonds_signature(m):
    # simple signature: sorted list of (a1_elem,a2_elem,order)
    sig = []
    for b in m.GetBonds():
        a1 = b.GetBeginAtom().GetSymbol()
        a2 = b.GetEndAtom().GetSymbol()
        order = int(b.GetBondTypeAsDouble())
        pair = tuple(sorted([a1,a2])) + (order,)
        sig.append(pair)
    return sorted(sig)

def has_double_bond_stereo(m):
    # Check if any double bond has stereochemistry assigned (E/Z)
    for b in m.GetBonds():
        if b.GetBondType().name == 'DOUBLE':
            if b.GetStereo() != Chem.BondStereo.STEREONONE:
                return True
    return False

def has_tetrahedral_chirality(m):
    for a in m.GetAtoms():
        if a.HasProp('_CIPCode') or a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            return True
    return False

def compare_files(pred_path, gt_path):
    pred = load_and_sanitize(pred_path)
    gt = load_and_sanitize(gt_path)
    out = {}
    out['pred_loaded'] = pred is not None
    out['gt_loaded'] = gt is not None
    if not pred or not gt:
        out.update({k: None for k in ['molecular_formula','molecular_bonds','double_bond_stereochemistry','tetrahedral_chirality']})
        return out
    out['molecular_formula'] = mol_formula(pred) == mol_formula(gt)
    out['molecular_bonds'] = bonds_signature(pred) == bonds_signature(gt)
    out['double_bond_stereochemistry'] = has_double_bond_stereo(pred) == has_double_bond_stereo(gt)
    out['tetrahedral_chirality'] = has_tetrahedral_chirality(pred) == has_tetrahedral_chirality(gt)
    # also include computed formulas for debugging
    out['pred_formula'] = mol_formula(pred)
    out['gt_formula'] = mol_formula(gt)
    return out

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: pb_rdkit_checks.py <predicted_sdf> <gt_sdf>')
        sys.exit(2)
    res = compare_files(sys.argv[1], sys.argv[2])
    for k,v in res.items():
        print(f'{k}: {v}')