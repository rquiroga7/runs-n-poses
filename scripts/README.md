# Docking Benchmark Scripts

Scripts to run the full docking benchmark on the Runs N' Poses single-ligand dataset
(1426 systems) and compare across docking programs.

## Layout

```
01_prepare_receptor_pdbqt*.py   — Obabel-based receptor prep (CIF→PDBQT)
02_prepare_ligand_pdbqt*.py      — Obabel-based ligand prep (SDF→PDBQT)
03_run_vina.py                   — AutoDock Vina runner
03_run_vinardo.py                — Vinardo (2vinardo-mar5) runner
04_analyze_vinardo.py            — Vinardo-specific analysis from obabel prep
04_analyze_vina.py / _32.py      — Vina-specific analysis (exhaustiveness 8/32)
04_analyze_docking.py            — Generic analysis (vina / vinardo / adgpu families)
05_plot_figure_with_*.py         — Figure 1E reproduction for each method
06_analyze_benchmark.py          — Speedup-vs-redocking scatter plots + summary
07_investigate_vinardock_failures.py — Failure investigation report
run_full_benchmark.py            — Orchestrator (receptor,ligand,dock,analyze steps)
run_full_benchmark_bg.sh         — Background launcher for run_full_benchmark.py

prepare_ligand_meeko.py          — Meeko-based ligand prep (SDF→PDBT + PDBQT)
prepare_receptor_meeko.py        — Meeko-based receptor prep (PDB→PDBT + PDBQT)
prepare_symmetry_*.py            — Obabel-based prep for symmetry-corrected inputs
symmetry_receptors_log.csv       — Log of symmetry generation

run_qvina_w.py                   — QuickVina-W runner
run_quickvina2.py                — QuickVina2 (scip-qvina) runner
run_vina_gpu21.py                — Vina-GPU 2.1 runner
run_vina_cuda.py                 — Vina-CUDA runner
run_vinardo_symmetry.py          — Vinardo runner for symmetry-corrected inputs
run_autodock_gpu.py              — AutoDock-GPU runner (uses Meeko-prepared inputs)
run_vina_symmetry.py             — Vina runner for symmetry-corrected inputs

*_symmetry.py                    — Runners for the symmetry-corrected pipeline
```

## Prerequisites

- `obabel-25-07` — OpenBabel for file conversion
- `vina`, `qvina-w`, `QuickVina2` — in PATH or at configured paths
- `vinardock-26-{mar,04}` — in `~/.local/bin/`
- `AutoDock-GPU`, `autogrid4` — at configured paths
- `Vina-CUDA_v1.1`, `AutoDock-Vina-GPU-2.1` — at configured paths
- `/home/rquiroga/github/Meeko` — Meeko for PDBT/PDBQT prep
- `/home/rquiroga/anaconda3/envs/runs_n_poses/bin/ost` — OpenStructure
- RDKit conda env at `/home/rquiroga/anaconda3/envs/RDKIT`
- PoseBusters venv at `/home/rquiroga/Datasets/Posebusters/venv`

## Input Data

Symmetry-corrected receptors and obabel-prepared ligands already exist:

```
/home/rquiroga/Datasets/runs-n-poses-datasets/
  symmetry_corrected/        — PDB files of symmetry-corrected complexes (2178)
  symmetry_receptors_pdbqt/  — Obabel-prepared PDBQT receptors
  symmetry_receptors_pdbt/   — Obabel-prepared PDBT receptors
  symmetry_ligands_pdbqt/    — Obabel-prepared PDBQT ligands
  symmetry_ligands_pdbt/     — Obabel-prepared PDBT ligands
  ground_truth/              — Original PLINDER ground truth (1426 single-ligand)
```

The benchmark runs on the 1426 single-ligand systems listed in
`scripts/single_ligand_systems_symmetry.txt`.

## Methods

| Method | Engine | Receptor Prep | Ligand Prep | Executable |
|--------|--------|---------------|-------------|------------|
| autodock_vina_8 | Vina 1.2 | obabel | obabel | `vina` (exh=8) |
| autodock_vina_32 | Vina 1.2 | obabel | obabel | `vina` (exh=32) |
| qvina_w | QuickVina-W | obabel | obabel | `qvina-w` |
| quickvina2 | QuickVina2 | obabel | obabel | `QuickVina2` |
| vina_gpu21 | Vina-GPU 2.1 | obabel | obabel | `AutoDock-Vina-GPU-2-1` |
| vina_cuda | Vina-CUDA 1.1 | obabel | obabel | `Vina-CUDA_v1.1` |
| autodock_gpu | AutoDock-GPU | Meeko | Meeko | `AutoDock-GPU` + `autogrid4` |
| vina_meeko | Vina 1.2 | Meeko | Meeko | `vina` (exh=8) |
| vinardock_2vinardo | 2vinardo-mar5 | obabel | obabel | `vinardock-26-mar` |
| vinardock_meeko | 2vinardo-mar5 | Meeko PDBT | Meeko PDBT | `vinardock-26-mar` |

## Running the Full Benchmark

```bash
# Background launcher (recommended):
bash scripts/run_full_benchmark_bg.sh

# Or run the orchestrator directly:
.venv/bin/python scripts/run_full_benchmark.py \
    --steps receptor,ligand,dock,analyze \
    --threads 8 --resume

# The orchestrator runs:
#   1. prepare_receptor_meeko.py   → meeko_receptors_{pdbqt,pdbt}/
#   2. prepare_ligand_meeko.py     → meeko_ligands_{pdbqt,pdbt}/
#   3. All 10 docking programs
#   4. 04_analyze_docking.py       → predictions/<method>.csv + posebusters/<method>.csv
```

Logs: `scripts/logs/benchmark/benchmark.log`, `scripts/logs/benchmark/pid`.

<rant>
Vina-GPU 2.1 and Vina-CUDA are known to fail on many systems (large boxes,
GPU resource limits). These failures are logged and skipped but expected.
</rant>

## Analyzing Results

```bash
# Speedup vs redocking% scatter + bar charts + failure analysis:
.venv/bin/python scripts/06_analyze_benchmark.py

# In-depth vinardock failure report:
.venv/bin/python scripts/07_investigate_vinardock_failures.py
```

Output goes to `figures/benchmark_*.png` and `figures/vinardock_failure_investigation.md`.

## Output Files

- `predictions/<method>.csv` — Prediction accuracy metrics
- `posebusters_results/<method>.csv` — PoseBusters physical plausibility
- `examples/analysis/<method>/` — OST comparison JSONs (cached)
- `figures/benchmark_*.png` — Benchmark plots
- `figures/vinardock_failure_*.csv|md` — Vinardock failure analysis
