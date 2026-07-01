#!/usr/bin/env python
"""
Step 2: Generate symmetry-corrected receptor structures for runs-n-poses.

Adapted from generate_sym_BETA.py for the runs-n-poses dataset.

Reads system.cif from ground_truth, identifies the proper (drug-like) ligand,
generates crystal symmetry mates, and outputs corrected receptor + ligand PDB files.

Usage:
    python prepare_symmetry_receptors.py \
        --system-info scripts/systems_for_symmetry_docking.csv \
        --ground-truth-dir runs-n-poses-datasets/ground_truth \
        --output-dir runs-n-poses-datasets/symmetry_corrected
"""

import os
import sys
import gemmi
import numpy as np
import traceback
import argparse
import csv
import gc
import urllib.request
import shutil

# --- Constants ---
EXCLUDED_RESIDUES = {'HOH', 'DOD', 'H2O', 'WAT'}
CONTACT_DISTANCE = 6.0
LIGAND_PROXIMITY = 8.0

AMINO_ACIDS = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY',
    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL', 'SEC', 'PYL', 'ACE', 'HCY',
    'CME', 'CYX', 'DPR', 'DTR', 'GLH', 'HID', 'HIE', 'HIM',
    'HSD', 'HSE', 'HSP', 'LYN', 'MSE', 'NLE', 'NVA', 'OCS',
    'ORN', 'PCA', 'PTR', 'SCY', 'SEP', 'TPO', 'TPQ', 'TRQ',
    'TYS', 'XLE', 'XLY', 'XME', 'XPE', 'XTR', 'XTY', 'XVA',
    'XYS', 'XYP',
}

# Ions that are common but should be removed from receptor if far from ligand
IONS_TO_REMOVE = {'NA', 'K', 'CL', 'BR', 'I', 'F'}

# RCSB download URL template
RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"


def get_pdb_crystal_params(pdb_id, cache_dir):
    """Get crystal parameters (cell, spacegroup) for a PDB ID.
    
    Downloads the PDB CIF from RCSB if not already cached, then extracts
    cell dimensions and space group.
    
    Returns (cell_obj, spacegroup_hm) or raises on failure.
    """
    pdb_lower = pdb_id.lower()
    # PDB CIFs are in subdirectories named by middle 2 chars
    subdir = pdb_lower[1:3]
    cif_path = os.path.join(cache_dir, subdir, f"{pdb_lower}.cif")

    if not os.path.exists(cif_path):
        os.makedirs(os.path.dirname(cif_path), exist_ok=True)
        url = RCSB_CIF_URL.format(pdb_id=pdb_lower)
        try:
            print(f"  Downloading {url}")
            urllib.request.urlretrieve(url, cif_path)
        except Exception as e:
            raise RuntimeError(f"Failed to download {url}: {e}")

    # Read just the metadata (cell + spacegroup) from the CIF without full parsing
    doc = gemmi.cif.read(cif_path)
    block = doc.sole_block()

    # Extract cell parameters (strip quotes from mmCIF string values)
    def _cif_float(tag):
        v = block.find_value(tag)
        if v is None:
            return 1.0 if "length" in tag else 90.0
        return float(str(v).strip().strip("'\""))

    cell_a = _cif_float("_cell.length_a")
    cell_b = _cif_float("_cell.length_b")
    cell_c = _cif_float("_cell.length_c")
    cell_alpha = _cif_float("_cell.angle_alpha")
    cell_beta = _cif_float("_cell.angle_beta")
    cell_gamma = _cif_float("_cell.angle_gamma")

    cell = gemmi.UnitCell(cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma)

    # Extract space group (try multiple mmCIF tags)
    sg_name = None
    for tag in ["_symmetry.space_group_name_H-M",
                "_symmetry_space_group_name_H-M",
                "_space_group.name_h-m_alt",
                "_symmetry.space_group_name_h-m",
                "_space_group_name_h-m_alt"]:
        val = block.find_value(tag)
        if val:
            sg_name = str(val).strip()
            break

    if not sg_name:
        raise RuntimeError(f"No space group found in {cif_path}")

    # Strip surrounding quotes (mmCIF convention for string values)
    sg_name = sg_name.strip().strip("'\"")

    return cell, sg_name


def is_ligand_residue(residue, ligand_asym_ids):
    """Check if a residue belongs to the proper ligand based on its subchain."""
    if not ligand_asym_ids:
        return False
    return residue.subchain in ligand_asym_ids


def is_protein_residue(residue):
    """Check if a residue is a standard or modified amino acid."""
    return residue.name in AMINO_ACIDS


def clear_anisou(structure):
    """Clear anisotropic B-factors from all atoms."""
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.aniso = gemmi.SMat33f(0, 0, 0, 0, 0, 0)


def assign_short_chain_ids(structure):
    """Rename all chains to single-letter PDB-compatible IDs (A-Z, then AA-ZZ).
    
    Returns the input structure with renamed chains (modified in-place).
    """
    # Collect existing chain names
    chains = list(structure[0])
    if not chains:
        return structure

    # Build mapping
    mapping = {}
    counter = 0
    for chain in chains:
        while True:
            if counter < 26:
                new_id = chr(65 + counter)
            else:
                first = chr(65 + (counter // 26) - 1)
                second = chr(65 + (counter % 26))
                new_id = f"{first}{second}"
            counter += 1
            if new_id not in mapping.values():
                break
        mapping[chain.name] = new_id

    # Apply mapping
    for chain in chains:
        chain.name = mapping[chain.name]
    
    return structure


def calculate_bounding_sphere(model):
    """Calculate the bounding sphere of a model (center + radius)."""
    x_sum, y_sum, z_sum = 0.0, 0.0, 0.0
    atom_count = 0

    for chain in model:
        for residue in chain:
            if residue.name in EXCLUDED_RESIDUES:
                continue
            for atom in residue:
                if atom.is_hydrogen():
                    continue
                x_sum += atom.pos.x
                y_sum += atom.pos.y
                z_sum += atom.pos.z
                atom_count += 1

    if atom_count == 0:
        return None

    center_x = x_sum / atom_count
    center_y = y_sum / atom_count
    center_z = z_sum / atom_count

    radius = 0.0
    for chain in model:
        for residue in chain:
            if residue.name in EXCLUDED_RESIDUES:
                continue
            for atom in residue:
                if atom.is_hydrogen():
                    continue
                dx = atom.pos.x - center_x
                dy = atom.pos.y - center_y
                dz = atom.pos.z - center_z
                dist = np.sqrt(dx*dx + dy*dy + dz*dz)
                radius = max(radius, dist)

    return (center_x, center_y, center_z, radius)


def can_spheres_intersect(sphere1, sphere2, max_distance):
    """Check if two bounding spheres could intersect within given distance."""
    dx = sphere1[0] - sphere2[0]
    dy = sphere1[1] - sphere2[1]
    dz = sphere1[2] - sphere2[2]
    center_dist = np.sqrt(dx*dx + dy*dy + dz*dz)
    return center_dist <= (sphere1[3] + sphere2[3] + max_distance)


def transform_sphere(sphere, sym_op, i, j, k, cell):
    """Transform a bounding sphere by symmetry operation and translation."""
    center_frac = cell.fractionalize(gemmi.Position(*sphere[0:3]))
    new_center = sym_op.apply_to_xyz([center_frac.x, center_frac.y, center_frac.z])
    new_center_orth = cell.orthogonalize(
        gemmi.Fractional(new_center[0] + i, new_center[1] + j, new_center[2] + k)
    )
    return (new_center_orth.x, new_center_orth.y, new_center_orth.z, sphere[3])


def has_chain_ligand_contacts(chain, ligand_st, cell, contact_distance):
    """Check if any non-water atom in chain contacts the ligand."""
    try:
        prot_positions = []
        for residue in chain:
            if residue.name in EXCLUDED_RESIDUES:
                continue
            for atom in residue:
                if atom.is_hydrogen():
                    continue
                prot_positions.append(atom.pos)

        if not prot_positions:
            return False

        for lig_chain in ligand_st[0]:
            for residue in lig_chain:
                for atom in residue:
                    if atom.is_hydrogen():
                        continue
                    for prot_pos in prot_positions:
                        if atom.pos.dist(prot_pos) <= contact_distance:
                            return True
        return False
    except Exception:
        return False


def apply_symmetry_to_chain(chain, sym_op, translation, cell, ligand_asym_ids=None):
    """Apply a symmetry operation and translation to a chain.
    
    Ligand residues get renamed to DpL; others keep original names.
    """
    new_chain = gemmi.Chain(chain.name)
    for residue in chain:
        new_residue = gemmi.Residue()
        if ligand_asym_ids and is_ligand_residue(residue, ligand_asym_ids):
            new_residue.name = "DpL"
        else:
            new_residue.name = residue.name
        new_residue.seqid = residue.seqid
        new_residue.subchain = residue.subchain
        new_residue.het_flag = residue.het_flag

        for atom in residue:
            if atom.is_hydrogen():
                continue
            new_atom = gemmi.Atom()
            new_atom.name = atom.name
            new_atom.element = atom.element
            if hasattr(atom, 'b_iso'):
                new_atom.b_iso = atom.b_iso
            if hasattr(atom, 'occ'):
                new_atom.occ = atom.occ

            pos = atom.pos
            frac = cell.fractionalize(pos)
            new_frac = sym_op.apply_to_xyz([frac.x, frac.y, frac.z])
            new_pos = cell.orthogonalize(gemmi.Fractional(
                new_frac[0] + translation[0],
                new_frac[1] + translation[1],
                new_frac[2] + translation[2]
            ))
            new_atom.pos = new_pos
            new_residue.add_atom(new_atom)

        new_chain.add_residue(new_residue)
    return new_chain


def build_crystal_environment(structure, ligand_asym_ids, cell, contact_distance, max_buffer_distance):
    """Build crystal environment: find symmetry operations with ligand contacts.
    
    Returns:
        candidate_operations: list of {"sym_op_idx": int, "translation": (i,j,k)}
        all_contact_distances: list of float
    """
    asu = structure[0]
    sg = gemmi.SpaceGroup(structure.spacegroup_hm)
    operations = list(sg.operations())

    asu_sphere = calculate_bounding_sphere(asu)
    if asu_sphere is None:
        return [], []

    # Add ligand atoms to neighbor search
    ns_asu = gemmi.NeighborSearch(asu, cell, contact_distance + max_buffer_distance)
    for chain_idx, chain in enumerate(asu):
        for res_idx, residue in enumerate(chain):
            if not is_ligand_residue(residue, ligand_asym_ids):
                continue
            for atom_idx, atom in enumerate(residue):
                if atom.is_hydrogen():
                    continue
                ns_asu.add_atom(atom, chain_idx, res_idx, atom_idx)

    # Determine search range
    a, b, c = cell.a, cell.b, cell.c
    search_radius = asu_sphere[3] + contact_distance + max_buffer_distance
    n_a = max(1, int(np.ceil(search_radius / a)))
    n_b = max(1, int(np.ceil(search_radius / b)))
    n_c = max(1, int(np.ceil(search_radius / c)))

    candidate_operations = []
    all_contact_distances = []

    for sym_idx, sym_op in enumerate(operations):
        if sym_idx == 0:
            continue  # skip identity
        for i in range(-n_a, n_a + 1):
            for j in range(-n_b, n_b + 1):
                for k in range(-n_c, n_c + 1):
                    new_sphere = transform_sphere(asu_sphere, sym_op, i, j, k, cell)
                    if not can_spheres_intersect(asu_sphere, new_sphere, contact_distance):
                        continue

                    has_contact = False
                    min_distance = float('inf')
                    contact_distances = []

                    for chain in asu:
                        for residue in chain:
                            if residue.name in EXCLUDED_RESIDUES:
                                continue
                            for atom in residue:
                                if atom.is_hydrogen():
                                    continue
                                pos = atom.pos
                                frac = cell.fractionalize(pos)
                                new_frac = sym_op.apply_to_xyz([frac.x, frac.y, frac.z])
                                new_pos = cell.orthogonalize(gemmi.Fractional(
                                    new_frac[0] + i, new_frac[1] + j, new_frac[2] + k
                                ))

                                results = ns_asu.find_atoms(new_pos, radius=contact_distance)
                                for mark in results:
                                    cra = mark.to_cra(asu)
                                    dist = new_pos.dist(cra.atom.pos)
                                    if dist <= contact_distance:
                                        contact_distances.append(dist)
                                    if dist < min_distance:
                                        min_distance = dist
                                        has_contact = dist <= contact_distance

                    if has_contact:
                        candidate_operations.append({
                            "sym_op_idx": sym_idx,
                            "translation": (i, j, k),
                            "min_distance": min_distance,
                        })
                        all_contact_distances.extend(contact_distances)

    candidate_operations.sort(key=lambda x: x["min_distance"])
    result = [{"sym_op_idx": op["sym_op_idx"], "translation": op["translation"]}
              for op in candidate_operations]
    return result, all_contact_distances


def read_system_info(csv_path):
    """Read the system info CSV: system_id -> proper_ligand_chain."""
    systems = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            systems.append((row["system_id"], row["proper_ligand_chain"]))
    return systems





def process_system(sys_id, proper_chain, ground_truth_dir, output_dir,
                   pdb_cache_dir, contact_distance, ligand_proximity):
    """Process a single system: generate symmetry-corrected receptor."""
    sys_out_dir = os.path.join(output_dir, sys_id)
    os.makedirs(sys_out_dir, exist_ok=True)

    ligand_pdb = os.path.join(sys_out_dir, f"{sys_id}_ligand.pdb")
    receptor_pdb = os.path.join(sys_out_dir, f"{sys_id}_receptor_symm.pdb")

    # Skip if both outputs already exist
    if os.path.exists(ligand_pdb) and os.path.exists(receptor_pdb):
        if os.path.getsize(ligand_pdb) > 0 and os.path.getsize(receptor_pdb) > 0:
            return "skip"

    system_cif = os.path.join(ground_truth_dir, sys_id, "system.cif")
    if not os.path.exists(system_cif):
        return f"no_system_cif: {system_cif}"

    print(f"\nProcessing {sys_id} (proper chain: {proper_chain})")

    # Extract PDB ID from system ID (first 4 characters)
    pdb_id = sys_id[:4]

    # Get crystal parameters from the original PDB CIF
    try:
        cell, sg_name = get_pdb_crystal_params(pdb_id, pdb_cache_dir)
    except Exception as e:
        return f"no_crystal_params: {e}"

    # 1. Read structure and override crystal parameters
    structure = gemmi.read_structure(system_cif)
    structure.cell = cell
    structure.spacegroup_hm = sg_name

    # 2. Build filtered structure: remove waters, distant non-proper ligands.
    # gemmi doesn't support in-place residue removal, so we build a new model.
    proper_positions = []
    for chain in structure[0]:
        for residue in chain:
            if is_ligand_residue(residue, {proper_chain}):
                for atom in residue:
                    if not atom.is_hydrogen():
                        proper_positions.append(atom.pos)

    filtered_model = gemmi.Model(0)
    for chain in structure[0]:
        new_chain = gemmi.Chain(chain.name)
        for residue in chain:
            # Skip waters
            if residue.name in EXCLUDED_RESIDUES:
                continue
            # Always keep proper ligand residues
            if is_ligand_residue(residue, {proper_chain}):
                new_chain.add_residue(residue)
                continue
            # Always keep protein
            if is_protein_residue(residue):
                new_chain.add_residue(residue)
                continue
            # Non-protein, non-proper-ligand: keep only if within proximity
            if proper_positions:
                min_dist = float('inf')
                for atom in residue:
                    if atom.is_hydrogen():
                        continue
                    for pp in proper_positions:
                        d = atom.pos.dist(pp)
                        if d < min_dist:
                            min_dist = d
                if min_dist <= ligand_proximity:
                    new_chain.add_residue(residue)
                # else: skip distant non-proper ligand
            else:
                new_chain.add_residue(residue)
        if len(new_chain):
            filtered_model.add_chain(new_chain)

    # Replace structure with filtered version
    structure = gemmi.Structure()
    structure.cell = cell
    structure.spacegroup_hm = sg_name
    structure.add_model(filtered_model)

    # 3. Split into proper ligand and receptor
    ligand_st = gemmi.Structure()
    ligand_st.cell = cell
    ligand_st.spacegroup_hm = structure.spacegroup_hm
    lig_model = gemmi.Model(0)

    receptor_st = gemmi.Structure()
    receptor_st.cell = cell
    receptor_st.spacegroup_hm = structure.spacegroup_hm
    rec_model = gemmi.Model(0)

    for chain in structure[0]:
        lig_chain = gemmi.Chain(chain.name)
        rec_chain = gemmi.Chain(chain.name)
        for residue in chain:
            if is_ligand_residue(residue, {proper_chain}):
                lig_chain.add_residue(residue)
            else:
                rec_chain.add_residue(residue)
        if len(lig_chain):
            lig_model.add_chain(lig_chain)
        if len(rec_chain):
            rec_model.add_chain(rec_chain)

    ligand_st.add_model(lig_model)
    receptor_st.add_model(rec_model)

    if len(ligand_st[0]) == 0:
        return "no_ligand_found"
    if len(receptor_st[0]) == 0:
        return "no_receptor_chains"

    # 4. Write proper ligand PDB (with short chain IDs)
    clear_anisou(ligand_st)
    assign_short_chain_ids(ligand_st)
    ligand_st.write_pdb(ligand_pdb)
    print(f"  Ligand written: {ligand_pdb}")

    # 5. Build crystal environment
    sg_name = structure.spacegroup_hm
    if not sg_name:
        return "no_spacegroup"

    ligand_asym_ids = {proper_chain}
    symm_ops, contact_dists = build_crystal_environment(
        structure, ligand_asym_ids, cell, contact_distance, max_buffer_distance=5.0
    )

    has_contacts = len(symm_ops) > 0
    if not has_contacts:
        print(f"  No symmetry contacts found — writing non-symm receptor only")
        clear_anisou(receptor_st)
        assign_short_chain_ids(receptor_st)
        receptor_st.write_pdb(receptor_pdb)
        return "no_symmetry_contacts"

    print(f"  Found {len(symm_ops)} symmetry operations with contacts")

    # 6. Generate symmetry mates
    sg_ops = list(gemmi.SpaceGroup(sg_name).operations())
    symm_input = structure.clone()
    crystal_model = gemmi.Model(0)
    used_chain_ids = {chain.name for chain in symm_input[0]}
    chain_id_counter = len(used_chain_ids)

    for op_info in symm_ops:
        sym_op = sg_ops[op_info["sym_op_idx"]]
        translation = op_info["translation"]
        orig_chain_names = [chain.name for chain in symm_input[0]]

        for orig_name in orig_chain_names:
            orig_chain = next(c for c in symm_input[0] if c.name == orig_name)

            # Generate unique chain ID
            while True:
                if chain_id_counter < 26:
                    new_chain_id = chr(65 + chain_id_counter)
                else:
                    first = chr(65 + (chain_id_counter // 26) - 1)
                    second = chr(65 + (chain_id_counter % 26))
                    new_chain_id = f"{first}{second}"
                if new_chain_id not in used_chain_ids:
                    break
                chain_id_counter += 1

            new_chain = apply_symmetry_to_chain(
                orig_chain, sym_op, translation, cell,
                ligand_asym_ids=ligand_asym_ids
            )
            if new_chain is None or len(new_chain) == 0:
                continue

            new_chain.name = new_chain_id

            # Check if original chain is a proper ligand chain
            orig_is_ligand = any(
                is_ligand_residue(r, ligand_asym_ids) for r in orig_chain
            )

            if orig_is_ligand:
                # Ligand symmetry mate: always add, residues already renamed to DpL
                crystal_model.add_chain(new_chain)
                used_chain_ids.add(new_chain_id)
                chain_id_counter += 1
            else:
                # Receptor chain: only add if it contacts the proper ligand
                if has_chain_ligand_contacts(new_chain, ligand_st, cell, contact_distance):
                    crystal_model.add_chain(new_chain)
                    used_chain_ids.add(new_chain_id)
                    chain_id_counter += 1

        gc.collect()

    # 7. Build final symmetry structure
    symm_st = gemmi.Structure()
    symm_st.cell = cell
    symm_st.spacegroup_hm = structure.spacegroup_hm
    symm_st.add_model(crystal_model)

    # Add original receptor chains
    for chain in receptor_st[0]:
        symm_st[0].add_chain(chain)

    # 8. Write receptor output (with short chain IDs)
    clear_anisou(symm_st)
    assign_short_chain_ids(symm_st)
    symm_st.write_pdb(receptor_pdb)
    print(f"  Receptor (symm) written: {receptor_pdb}")
    print(f"  Total chains in output: {len(symm_st[0])}")

    # Cleanup
    del structure, ligand_st, receptor_st, symm_st, symm_input, crystal_model
    gc.collect()

    return "ok"


def main():
    parser = argparse.ArgumentParser(
        description="Generate symmetry-corrected receptors for runs-n-poses"
    )
    parser.add_argument(
        "--system-info",
        default=os.path.join(os.path.dirname(__file__), "systems_for_symmetry_docking.csv"),
        help="CSV with system_id,proper_ligand_chain"
    )
    parser.add_argument(
        "--ground-truth-dir",
        default=os.path.join(os.path.dirname(__file__), "..",
                            "runs-n-poses-datasets", "ground_truth"),
        help="Path to ground_truth directory"
    )
    parser.add_argument(
        "--output-dir",
        default=os.path.join(os.path.dirname(__file__), "..",
                            "runs-n-poses-datasets", "symmetry_corrected"),
        help="Output directory for symmetry-corrected PDB files"
    )
    parser.add_argument(
        "--pdb-cache",
        default=os.path.join(os.path.dirname(__file__), "..",
                            "runs-n-poses-datasets", "pdb_cache"),
        help="Directory to cache downloaded PDB CIF files"
    )
    parser.add_argument(
        "--contact-distance", type=float, default=CONTACT_DISTANCE,
        help=f"Contact distance for symmetry detection (default: {CONTACT_DISTANCE})"
    )
    parser.add_argument(
        "--ligand-proximity", type=float, default=LIGAND_PROXIMITY,
        help=f"Max distance for keeping non-proper ligands (default: {LIGAND_PROXIMITY})"
    )
    parser.add_argument(
        "--system-id", type=str, default=None,
        help="Process only a specific system ID (for testing)"
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="Skip systems where output files already exist"
    )
    args = parser.parse_args()

    ground_truth_dir = os.path.abspath(args.ground_truth_dir)
    output_dir = os.path.abspath(args.output_dir)
    pdb_cache_dir = os.path.abspath(args.pdb_cache)
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(pdb_cache_dir, exist_ok=True)

    systems = read_system_info(args.system_info)
    if args.system_id:
        systems = [s for s in systems if s[0] == args.system_id]

    print(f"Processing {len(systems)} systems")
    print(f"  Contact distance: {args.contact_distance} Å")
    print(f"  Ligand proximity: {args.ligand_proximity} Å")
    print(f"  PDB cache: {pdb_cache_dir}")
    print(f"  Output: {output_dir}")

    # Track results
    results = {"ok": 0, "skip": 0, "no_symmetry_contacts": 0, "failed": 0}
    log_entries = []

    for sys_id, proper_chain in systems:
        try:
            status = process_system(
                sys_id, proper_chain, ground_truth_dir, output_dir,
                pdb_cache_dir, args.contact_distance, args.ligand_proximity
            )
            if status == "ok":
                results["ok"] += 1
            elif status == "skip":
                results["skip"] += 1
                continue
            elif status.startswith("no_symmetry"):
                results["no_symmetry_contacts"] += 1
                log_entries.append((sys_id, status))
            else:
                results["failed"] += 1
                log_entries.append((sys_id, status))
            print(f"  -> {status}")
        except Exception as e:
            trace = traceback.format_exc()
            results["failed"] += 1
            log_entries.append((sys_id, f"{type(e).__name__}: {e}"))
            print(f"  ERROR: {e}")
            print(trace)
        gc.collect()

    print(f"\n{'='*60}")
    print(f"Results: OK={results['ok']}, Skip={results['skip']}, "
          f"NoSymmContacts={results['no_symmetry_contacts']}, Failed={results['failed']}")
    print(f"{'='*60}")

    # Write log
    log_path = os.path.join(os.path.dirname(__file__), "symmetry_receptors_log.csv")
    with open(log_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["system_id", "status"])
        for sys_id, status in log_entries:
            writer.writerow([sys_id, status])
    print(f"Log written to {log_path}")


if __name__ == "__main__":
    main()
