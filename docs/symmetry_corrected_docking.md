# Symmetry-Corrected Docking for Runs-N-Poses

## Objective

Generate biologically relevant symmetry-corrected receptor structures for the
runs-n-poses dataset and re-dock the proper (drug-like) ligands using AutoDock
Vina (exhaustiveness=32) and VinardoCK (2vinardo scoring).

This accounts for crystal packing effects where symmetry-related protein chains
form part of the ligand binding site.

---

## Dataset

| Category | Count |
|---|---|
| Single-ligand systems | 1426 |
| Multi-ligand systems with exactly 1 proper ligand | 752 |
| **Total target systems** | **2178** |

Multi-ligand systems included only when one ligand is drug-like
(`ligand_is_proper=True`) and all other ligands are non-drug (metals, ions,
co-solvents such as ZN, MG, SO4, GOL, EDO, DMS, ACT, etc.).

---

## Pipeline

### Step 1 — Generate system list

**Script**: `scripts/generate_system_list.py`

**Input**: `runs-n-poses-datasets/annotations.csv`

**Output files**:
- `scripts/systems_for_symmetry_docking.txt` — one system_id per line (2178 entries)
- `scripts/systems_for_symmetry_docking.csv` — `system_id,proper_ligand_chain`

**Logic**:
- Single-ligand systems (`num_ligand_chains == 1`): the only chain is proper
- Multi-ligand systems with exactly 1 proper (`num_ligand_chains > 1` and
  `sum(ligand_is_proper) == 1`): identify proper chain from `ligand_instance_chain`

---

### Step 2 — Generate symmetry-corrected receptors

**Script**: `scripts/prepare_symmetry_receptors.py`

**Parameters**:

| Parameter | Value |
|---|---|
| Contact distance (symmetry mates) | 6.0 Å |
| Ligand proximity (keep non-proper ligands) | 8.0 Å |
| Waters | Removed (HOH, DOD, H2O, WAT) |
| Non-proper ligands >8Å from proper ligand | Removed |
| Symmetry ligand copies | Renamed to DpL (prevents multiple binding sites) |

**Per-system workflow**:

1. Read `ground_truth/<sys_id>/system.cif`
2. Identify proper ligand (by chain from CSV), confirm via SDF coordinate match
3. Keep non-proper ligands within 8Å of proper ligand atoms → merge into receptor
4. Build crystal environment: enumerate space group ops + unit cell translations;
   find those where the transformed ASU has atoms within 6Å of the proper ligand
5. Generate symmetry mates:
   - All ligand atoms in symmetry copies → rename to DpL
   - Only add protein symmetry mates that directly contact the proper ligand
   - Always add non-protein symmetry mates (they become DpL)
6. Cleanup: remove waters
7. Output files:

```
symmetry_corrected/<sys_id>/
  <sys_id>_ligand.pdb              # The proper drug-like ligand
  <sys_id>_receptor_symm.pdb       # Receptor + symmetry mates (docking target)
```

---

### Step 3 — Prepare docking input files

Three scripts, each loops over the 2178 systems.

#### 3a. Receptor PDBT for Vinardo

**Script**: `scripts/prepare_receptor_pdbt.py`

```
symmetry_corrected/<sys_id>/<sys_id>_receptor_symm.pdb
  → obabel -opdbt -h -xc -xr
  → symmetry_receptors_pdbt/<sys_id>/<sys_id>_receptor.pdbt
```

#### 3b. Receptor PDBQT for Vina

**Script**: `scripts/prepare_receptor_pdbqt.py`

Same as 3a but `-opdbqt`. If the output contains `ROOT`/`ENDROOT` tags (which Vina
may reject), a round-trip PDBQT→PDB→PDBQT cleanup is performed automatically.

```
symmetry_corrected/<sys_id>/<sys_id>_receptor_symm.pdb
  → obabel -opdbqt -h -xc -xr  (+ ROOT cleanup if needed)
  → symmetry_receptors_pdbqt/<sys_id>/<sys_id>_receptor.pdbqt
```

#### 3c. Ligands for both methods

**Script**: `scripts/prepare_symmetry_ligands.py`

```
ground_truth/<sys_id>/ligand_files/<chain>.sdf
  → obabel -opdbt -h              → symmetry_ligands_pdbt/<sys_id>/<chain>.pdbt
  → obabel -opdbqt -h             → symmetry_ligands_pdbqt/<sys_id>/<chain>.pdbqt
```

---

### Step 4 — Run docking

#### 4a. AutoDock Vina (exhaustiveness=32)

**Script**: `scripts/run_vina_symmetry.py`

- Receptors: `symmetry_receptors_pdbqt/`
- Ligands: `symmetry_ligands_pdbqt/`
- Exhaustiveness: 32
- Output: `autodock_vina_32_symmetry/<sys_id>/<sys_id>_<chain>/out.pdbqt`
- Autobox from receptor coordinates + 10Å margin

#### 4b. VinardoCK (2vinardo scoring)

**Script**: `scripts/run_vinardo_symmetry.py`

- Receptors: `symmetry_receptors_pdbt/`
- Ligands: `symmetry_ligands_pdbt/`
- Scoring: `--scoring 2vinardo`
- Output: `vinardo_outputs_symmetry/<sys_id>/<sys_id>_<chain>/log.csv`
- Autobox from receptor coordinates

---

## Directory layout

```
runs-n-poses-datasets/
├── symmetry_corrected/                 # Step 2
│   └── <sys_id>/
│       ├── <sys_id>_ligand.pdb
│       └── <sys_id>_receptor_symm.pdb
├── symmetry_receptors_pdbt/            # Step 3a
│   └── <sys_id>/<sys_id>_receptor.pdbt
├── symmetry_receptors_pdbqt/           # Step 3b
│   └── <sys_id>/<sys_id>_receptor.pdbqt
├── symmetry_ligands_pdbt/              # Step 3c
│   └── <sys_id>/<chain>.pdbt
├── symmetry_ligands_pdbqt/             # Step 3c
│   └── <sys_id>/<chain>.pdbqt
├── autodock_vina_32_symmetry/          # Step 4a
│   └── <sys_id>/<sys_id>_<chain>/
│       ├── out.pdbqt
│       └── log.txt
└── vinardo_outputs_symmetry/           # Step 4b
    └── <sys_id>/<sys_id>_<chain>/
        ├── log.csv
        └── <chain>_dock.sdf
```

---

## Execution order

1. `generate_system_list.py` — create the 2178-entry system list
2. `prepare_symmetry_receptors.py` — symmetry generation for all systems
3a. `prepare_receptor_pdbt.py` — receptor PDBT conversion
3b. `prepare_receptor_pdbqt.py` — receptor PDBQT conversion (+ ROOT cleanup)
3c. `prepare_symmetry_ligands.py` — ligand PDBT + PDBQT for all systems
4a. `run_vina_symmetry.py` — Vina docking (exh=32)
4b. `run_vinardo_symmetry.py` — Vinardo docking

---

## Key design decisions

1. **Multi-ligand systems**: only include if exactly 1 proper (drug-like) ligand.
   Non-proper ligands (metals, ions, co-solvents) within 8Å of the proper ligand
   are kept in the receptor; farther ones are removed.
2. **Single binding site**: all ligand atoms in symmetry copies are renamed to DpL,
   ensuring only the original binding site is recognized by the docking program.
3. **Contact-based filtering**: only protein symmetry mates that have contacts
   within 6Å of the proper ligand are included, reducing file size and noise.
4. **Ligands from ground truth SDFs**: always use original crystal coordinates;
   never generate from SMILES.
5. **Receptor conversion**: obabel with `-h -xc -xr` handles hydrogen addition,
   charge removal, and problematic residue filtering.
6. **ROOT-tag cleanup**: PDBQT files with `ROOT`/`ENDROOT` tags are cleaned via
   obabel round-trip (PDBQT→PDB→PDBQT).

---

## Known risks & mitigations

| Risk | Mitigation |
|---|---|
| gemmi memory leaks on large CIFs | Explicit `del` + `gc.collect()` per iteration |
| Some systems may have no symmetry contacts | Log to tracking file; skip symmetry for those |
| obabel crashes on non-standard residues | `-xr` flag removes problem residues |
| Vinardo crashes on unsupported metal types | `clean_receptor_pdbt.py` post-processing if needed |
| Multi-ligand where non-proper IS the binding site partner | Rare; shows as poor docking scores downstream |

---

## Scripts

| Step | Script | Input | Output | Description |
|---|---|---|---|---|---|
| 1 | `generate_system_list.py` | `annotations.csv` | `systems_for_symmetry_docking.{txt,csv}` | Create 2178-entry system list |
| 2 | `prepare_symmetry_receptors.py` | `ground_truth/*/system.cif` + PDB CIF cache | `symmetry_corrected/<sys_id>/` | Generate symmetry-corrected receptors |
| 3a | `prepare_receptor_pdbt.py` | `symmetry_corrected/` | `symmetry_receptors_pdbt/` | Symm PDB → PDBT |
| 3b | `prepare_receptor_pdbqt.py` | `symmetry_corrected/` | `symmetry_receptors_pdbqt/` | Symm PDB → PDBQT (+ ROOT cleanup) |
| 3c | `prepare_symmetry_ligands.py` | `ground_truth/*/ligand_files/*.sdf` | `symmetry_ligands_{pdbt,pdbqt}/` | Ligand PDBT + PDBQT |
| 4a | `run_vina_symmetry.py` | receptors + ligands PDBQT | `autodock_vina_32_symmetry/` | Vina docking (exh=32) |
| 4b | `run_vinardo_symmetry.py` | receptors PDBT + ligands PDBT | `vinardo_outputs_symmetry/` | Vinardo docking (`--scoring 2vinardo`) |
| 4c | `run_qvina_w_symmetry.py` | receptors + ligands PDBQT | `qvina_w_symmetry/` | QVina-W v1.2 |
| 4d | `run_quickvina2_symmetry.py` | receptors + ligands PDBQT | `quickvina2_symmetry/` | QuickVina2 (scip-qvina) |
| 4e | `run_psovina2_symmetry.py` | receptors + ligands PDBQT | `psovina2_symmetry/` | PSOVina2 |
| 4f | `run_vinagpu21_symmetry.py` | receptors + ligands PDBQT | `vinagpu21_symmetry/` | Vina-GPU-2.1 |
| 4g | `run_vina_cuda_symmetry.py` | receptors + ligands PDBQT | `vina_cuda_symmetry/` | Vina-CUDA v1.1 |
| 4h | AutoDock-GPU (Scripps) | grid maps `.fld` | `adgpu_symmetry/` | Needs autogrid4 maps pre-generation |

## Pipeline status

| Step | Systems | Success | Fail | Notes |
|---|---|---|---|---|
| 1 — System list | 2178 | 2178 | 0 | 1426 single + 752 multi-ligand |
| 2 — Symmetry receptors | 2178 | 2177 | 1 | 8c3i__1__1.B__1.D (corrupted system.cif) |
| 3a — Receptor PDBT | 2178 | 2177 | 1 | 8c3i__1__1.B__1.D missing |
| 3b — Receptor PDBQT | 2178 | 2177 | 1 | 8c3i__1__1.B__1.D missing |
| 3c — Ligands | 2178 | 2178 | 0 | |
| 4a — Vina (exh=32) | 2178 | 1898 | 280 | 264 not started + 15 no output + 1 no receptor (8c3i) |
| 4b — Vinardo | 2178 | 2177 | 1 | 8c3i__1__1.B__1.D (corrupted system.cif) |
| 4c — QVina-W v1.2 | 2178 | Pending | — | |
| 4d — QuickVina2 (scip-qvina) | 2178 | Pending | — | |
| 4e — PSOVina2 | 2178 | Pending | — | |
| 4f — Vina-GPU-2.1 | 2178 | Pending | — | |
| 4g — Vina-CUDA v1.1 | 2178 | Pending | — | |
| 4h — AutoDock-GPU (Scripps) | 2178 | Pending | — | Needs autogrid4 maps pre-generation |
