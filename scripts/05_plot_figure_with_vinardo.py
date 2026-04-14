#!/usr/bin/env python
"""
Step 5: Reproduce Figure 1 Panel E from the paper with vinardock results added.

This script follows the exact same pattern as figures.ipynb:
1. Load annotations, predictions, and posebusters results
2. Process data with pivot_df (same as notebook)
3. Create common_subset_dfs_all["top"]
4. Call plotting.make_main_figure() with vinardock included

Usage:
    python 05_plot_figure_with_vinardo.py \
        --data-dir /path/to/data \
        --output /path/to/output_figure.png
"""

import os
os.environ["MPLBACKEND"] = "Agg"

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Add parent directory to path to import plotting module
sys.path.insert(0, str(Path(__file__).parent.parent))
import plotting


def load_data(data_dir: Path, methods: list):
    """
    Load all data following the exact same pattern as figures.ipynb.
    """
    # Load annotations
    annotated_df = pd.read_csv(data_dir / "annotations.csv")
    annotated_df["release_date"] = pd.to_datetime(annotated_df["release_date"])

    # Load posebusters results
    bust_dfs = {}
    for m in methods:
        filename = data_dir / "posebusters_results" / f"{m}.csv"
        if filename.exists():
            bust_df = pd.read_csv(filename)
            # Rename ligand_chain -> ligand_instance_chain if needed
            if "ligand_chain" in bust_df.columns and "ligand_instance_chain" not in bust_df.columns:
                bust_df["ligand_instance_chain"] = bust_df["ligand_chain"]
            # Compute pb_success from individual checks if not present
            if "pb_success" not in bust_df.columns:
                checks = ["sanitization", "inchi_convertible", "all_atoms_connected",
                          "bond_lengths", "bond_angles", "internal_steric_clash",
                          "aromatic_ring_flatness", "double_bond_flatness",
                          "internal_energy", "protein-ligand_maximum_distance",
                          "minimum_distance_to_protein"]
                available = [c for c in checks if c in bust_df.columns]
                if available:
                    bust_df["pb_success"] = bust_df[available].all(axis=1).astype(float)
                else:
                    bust_df["pb_success"] = -1
            bust_dfs[m] = bust_df

    # Load predictions
    full_datasets = {}
    for method in methods:
        filename = data_dir / "predictions" / f"{method}.csv"
        if not filename.exists():
            print(f"Warning: {filename} not found, skipping {method}")
            continue

        df = pd.read_csv(filename, low_memory=False)
        keep_columns = [
            "target",
            "ligand_instance_chain",
            "lddt_pli",
            "rmsd",
            "lddt_lp",
            "bb_rmsd",
            "seed",
            "sample",
            "ranking_score",
            "ligand_is_proper",
            "prot_lig_chain_iptm_average_lddt_pli",
            "prot_lig_chain_iptm_min_lddt_pli",
            "prot_lig_chain_iptm_max_lddt_pli",
            "lig_prot_chain_iptm_average_lddt_pli",
            "lig_prot_chain_iptm_min_lddt_pli",
            "lig_prot_chain_iptm_max_lddt_pli",
            "prot_lig_chain_iptm_average_rmsd",
            "prot_lig_chain_iptm_min_rmsd",
            "prot_lig_chain_iptm_max_rmsd",
            "lig_prot_chain_iptm_average_rmsd",
            "lig_prot_chain_iptm_min_rmsd",
            "lig_prot_chain_iptm_max_rmsd",
            "model_ligand_chain_lddt_pli",
            "model_ligand_chain_rmsd",
            "ligand_ccd_code",
            "model_ligand_smiles",
            "pred_pocket_f1",
        ]
        if "seed" not in df.columns:
            df["seed"] = 1
        if "sample" not in df.columns:
            df["sample"] = 1
        if "ranking_score" not in df.columns:
            df["ranking_score"] = 1
        if "lig_prot_chain_iptm_average_rmsd" not in df.columns:
            df["lig_prot_chain_iptm_average_rmsd"] = 1
        if "prot_lig_chain_iptm_average_rmsd" not in df.columns:
            df["prot_lig_chain_iptm_average_rmsd"] = 1
        if "pred_pocket_f1" not in df.columns:
            df["pred_pocket_f1"] = 1
        keep_columns = [c for c in keep_columns if c in df.columns]
        full_datasets[method] = (
            df[keep_columns].rename(columns={"target": "system_id"}).reset_index(drop=True)
        )
        full_datasets[method]["group_key"] = (
            full_datasets[method]["system_id"]
            + "__"
            + full_datasets[method]["ligand_instance_chain"]
        )
        full_datasets[method]["method"] = method
        full_datasets[method] = (
            full_datasets[method]
            .sort_values(by=["lddt_pli", "rmsd"], ascending=[False, True])
            .groupby(["group_key", "seed", "sample"])
            .head(1)
            .reset_index(drop=True)
        )
        # Merge PoseBusters results if available. If not, try to fallback
        # to the `pb_success` column inside the predictions CSV itself.
        if method in bust_dfs:
            full_datasets[method] = full_datasets[method].merge(
                bust_dfs[method][["system_id", "ligand_instance_chain", "pb_success"]],
                on=["system_id", "ligand_instance_chain"],
                how="left",
            )
            full_datasets[method]["pb_success"] = (
                full_datasets[method]["pb_success"].fillna(False).astype(float)
            )
        else:
            # Try to read pb_success from the predictions file we already loaded
            if "pb_success" in df.columns:
                bust_df = (
                    df[["target", "ligand_instance_chain", "pb_success"]]
                    .rename(columns={"target": "system_id"})
                )
                full_datasets[method] = full_datasets[method].merge(
                    bust_df[["system_id", "ligand_instance_chain", "pb_success"]],
                    on=["system_id", "ligand_instance_chain"],
                    how="left",
                )
                full_datasets[method]["pb_success"] = (
                    full_datasets[method]["pb_success"].fillna(False).astype(float)
                )
            else:
                full_datasets[method]["pb_success"] = -1

    return full_datasets, annotated_df


def pivot_df(df, annotated_df):
    """Exact copy from figures.ipynb."""
    df = df.pivot(
        index=[
            "group_key",
            "system_id",
            "ligand_is_proper",
            "ligand_instance_chain",
        ],
        columns="method",
        values=[
            "lddt_pli",
            "rmsd",
            "lddt_lp",
            "bb_rmsd",
            "pb_success",
            "pred_pocket_f1",
        ],
    ).reset_index()
    df.columns = [f"{col[0]}_{col[1]}" if len(col[1]) else col[0] for col in df.columns]
    df = df[df["ligand_is_proper"].fillna(False)].reset_index(drop=True)
    merge_columns = [col for col in annotated_df.columns if col not in df.columns]
    df = df.merge(
        annotated_df[["group_key"] + merge_columns], on="group_key", how="left"
    )
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Reproduce Figure 1 with vinardock results (follows figures.ipynb)"
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets",
        help="Path to data directory (contains annotations.csv, predictions/, posebusters_results/)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="/home/rquiroga/github/runs-n-poses/figures/fig_1_with_vinardo.png",
        help="Path to output figure file",
    )
    parser.add_argument(
        "--vinardo-method-name",
        type=str,
        default="vinardock_2vinardo",
        help="Method name for vinardo in predictions",
    )
    parser.add_argument(
        "--common-subset-only",
        action="store_true",
        help="Only include systems where ALL methods (including vinardo) have results",
    )

    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Replace Protenix with VinardoDock in the common subset
    base_methods = [m for m in plotting.COMMON_SUBSET_METHODS if m != "protenix"]
    base_methods.append(args.vinardo_method_name)

    # Always include all methods in the common subset
    all_methods = base_methods
    common_subset_methods = base_methods

    print(f"Loading data for methods: {all_methods}")
    print(f"Common subset methods: {common_subset_methods}")

    # Load data (exact same pattern as figures.ipynb)
    full_datasets, annotated_df = load_data(data_dir, all_methods)

    if not full_datasets:
        print("Error: No prediction data loaded")
        sys.exit(1)

    print(f"Loaded {len(full_datasets)} methods: {list(full_datasets.keys())}")

    # Process data (exact same pattern as figures.ipynb cell 9)
    top_dfs = {}
    rank_by = "lig_prot_chain_iptm_average_rmsd"
    
    for m in full_datasets:
        top_dfs[m] = (
            full_datasets[m]
            .sort_values(by=rank_by, ascending=False)
            .groupby(["system_id", "ligand_instance_chain"])
            .head(1)
        )

    results_df_top = pivot_df(pd.concat(top_dfs.values()), annotated_df)

    # Add best/worst/average metrics (same as notebook)
    # Extract method names from pivoted column names like "lddt_pli_af3"
    prefix = "lddt_pli_"
    methods_in_df = list(set([
        col[len(prefix):]
        for col in results_df_top.columns
        if col.startswith(prefix) and col[len(prefix):] not in ("max", "average")
    ]))
    # Sort to match the order in plotting.METHODS + vinardo
    methods_in_df = [m for m in base_methods if m in methods_in_df]

    results_df_top["lddt_pli_max"] = np.nanmax(
        results_df_top[
            [
                f"lddt_pli_{m}"
                for m in plotting.METHODS
                if f"lddt_pli_{m}" in results_df_top.columns
            ]
        ],
        axis=1,
    )
    results_df_top["rmsd_min"] = np.nanmin(
        results_df_top[
            [
                f"rmsd_{m}"
                for m in plotting.METHODS
                if f"rmsd_{m}" in results_df_top.columns
            ]
        ],
        axis=1,
    )
    results_df_top["lddt_pli_average"] = np.nanmedian(
        results_df_top[
            [
                f"lddt_pli_{m}"
                for m in plotting.METHODS
                if f"lddt_pli_{m}" in results_df_top.columns
            ]
        ],
        axis=1,
    )
    results_df_top["rmsd_average"] = np.nanmedian(
        results_df_top[
            [
                f"rmsd_{m}"
                for m in plotting.METHODS
                if f"rmsd_{m}" in results_df_top.columns
            ]
        ],
        axis=1,
    )

    # Create common subset (same as notebook)
    common_subset_df_top = (
        results_df_top
        .dropna(
            subset=[f"lddt_pli_{method}" for method in common_subset_methods]
            + ["sucos_shape"]
        )
        .reset_index(drop=True)
    )

    print(f"\nData summary:")
    print(f"  Total systems (top): {results_df_top['system_id'].nunique()}")
    print(f"  Common subset systems: {common_subset_df_top['system_id'].nunique()}")
    print(f"  Methods in results: {methods_in_df}")

    # Add vinardock color and shape to plotting module
    plotting.COLORS[args.vinardo_method_name] = "#FF69B4"  # Hot pink
    plotting.SHAPES[args.vinardo_method_name] = "^"  # Triangle
    plotting.NAME_MAPPING[args.vinardo_method_name] = "VinardoDock"

    # Create figure using plotting.make_main_figure (exact same call as figures.ipynb)
    print("\nCreating figure using plotting.make_main_figure()...")
    
    plotting.make_main_figure(
        common_subset_df_top,
        str(output_path),
        methods=methods_in_df,
    )

    print(f"\nFigure saved to: {output_path}")
    print("Done!")


if __name__ == "__main__":
    main()
