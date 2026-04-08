# Vinardo Docking Pipeline Scripts

This directory contains scripts to run the 2vinardo-mar5 docking pipeline and analyze the results in the context of the Runs N' Poses benchmark.

## Overview

The pipeline consists of 5 steps:

1. **Prepare receptor PDBQT files** - Convert receptor.cif files to PDBQT format
2. **Prepare ligand PDBQT files** - Convert ligand SMILES/SDF files to PDBQT format  
3. **Run 2vinardo-mar5 docking** - Execute the docking for each receptor-ligand pair
4. **Analyze results** - Score predictions and output in the same format as other methods
5. **Plot figures** - Create comparison plots including vinardock results

## Prerequisites

- `obabel-25-07` - OpenBabel for file format conversion
- `2vinardo-mar5_autobox` - The docking executable (must be in PATH)
- `ost` - OpenStructure for structure comparison (must be in PATH)
- Python 3.8+ with packages from `environment.yaml`

## Usage

### Step 1: Prepare Receptor PDBQT Files

```bash
python 01_prepare_receptor_pdbqt.py \
    --ground-truth-dir /home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth \
    --output-dir /home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_inputs/receptors \
    --annotations /home/rquiroga/Datasets/runs-n-poses-datasets/annotations.csv
```

This converts receptor.cif files to PDBQT format using `obabel-25-07` with `-xc -xr` flags to handle cofactors and metals properly.

### Step 2: Prepare Ligand PDBQT Files

```bash
python 02_prepare_ligand_pdbqt.py \
    --annotations /home/rquiroga/Datasets/runs-n-poses-datasets/annotations.csv \
    --ground-truth-dir /home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth \
    --output-dir /home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_inputs/ligands
```

This creates PDBQT files for each ligand from the SMILES strings in annotations.csv (or from ground truth SDF files if SMILES conversion fails).

### Step 3: Run 2vinardo-mar5 Docking

```bash
python 03_run_vinardo.py \
    --receptor-dir /home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_inputs/receptors \
    --ligand-dir /home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_inputs/ligands \
    --output-dir /home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_outputs \
    --config /home/rquiroga/github/runs-n-poses/config.fijo
```

This runs `2vinardo-mar5_autobox` for each receptor-ligand pair with autoboxing enabled.

### Step 4: Analyze Vinardo Results

```bash
python 04_analyze_vinardo.py \
    --vinardo-outputs /home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_outputs \
    --ground-truth /home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth \
    --annotations /home/rquiroga/Datasets/runs-n-poses-datasets/annotations.csv \
    --inputs-json /home/rquiroga/Datasets/runs-n-poses-datasets/inputs.json \
    --output /home/rquiroga/Datasets/runs-n-poses-datasets/predictions/vinardock_2vinardo.csv
```

This scores the vinardo predictions against ground truth using `ost compare-ligand-structures` (following the pattern in `examples/utils/analyze_models_*.sh` and `extract_scores.ipynb`) and outputs a CSV in the same format as the other prediction methods.

### Step 5: Plot Figures with Vinardo Results

```bash
python 05_plot_figure_with_vinardo.py \
    --predictions-dir /home/rquiroga/Datasets/runs-n-poses-datasets/predictions \
    --annotations /home/rquiroga/Datasets/runs-n-poses-datasets/annotations.csv \
    --output /home/rquiroga/github/runs-n-poses/figures/figure_1e_with_vinardo.png \
    --comparison-output /home/rquiroga/github/runs-n-poses/figures/figure_1e_comparison.png
```

This reproduces Figure 1 Panel E from the paper with vinardock results added.

## Running the Full Pipeline

To run all steps sequentially:

```bash
# Step 1
python 01_prepare_receptor_pdbqt.py

# Step 2
python 02_prepare_ligand_pdbqt.py

# Step 3
python 03_run_vinardo.py

# Step 4
python 04_analyze_vinardo.py

# Step 5
python 05_plot_figure_with_vinardo.py
```

Each script supports a `--system-id` flag to process only a specific system for testing:

```bash
python 01_prepare_receptor_pdbqt.py --system-id "5s9l__1__1.A__1.H_1.I"
python 02_prepare_ligand_pdbqt.py --system-id "5s9l__1__1.A__1.H_1.I"
python 03_run_vinardo.py --system-id "5s9l__1__1.A__1.H_1.I"
python 04_analyze_vinardo.py --system-id "5s9l__1__1.A__1.H_1.I"
```

## Output Files

- `vinardo_inputs/receptors/` - Receptor PDBQT files
- `vinardo_inputs/ligands/` - Ligand PDBQT files
- `vinardo_outputs/` - Raw docking output files
- `vinardo_analysis/` - Intermediate JSON scoring files
- `predictions/vinardock_2vinardo.csv` - Final predictions CSV
- `figures/figure_1e_with_vinardo.png` - Figure with vinardo results
- `figures/figure_1e_comparison.png` - Line plot comparison

## Notes

- The pipeline uses the same scoring methodology as the original Runs N' Poses benchmark
- Vinardo doesn't produce iPTM scores (these are set to None in the output)
- The `config.fijo` file contains the Vinardo scoring parameters and enables autobox mode
- All scripts backup existing files before overwriting (`.bkp` extension)
