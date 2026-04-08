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
            bust_dfs[m] = pd.read_csv(filename)

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
    parser.add_argument(
        "--no-vinardo-in-common-subset",
        action="store_true",
        help="Use original common subset methods but add vinardo if available (recommended)",
    )

    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Define methods - following plotting.COMMON_SUBSET_METHODS + vinardo
    base_methods = list(plotting.COMMON_SUBSET_METHODS)  # ["af3", "protenix", "chai", "boltz"]
    
    if args.no_vinardo_in_common_subset:
        # Common subset uses original methods, but vinardo is added for visualization if available
        all_methods = base_methods + [args.vinardo_method_name]
        common_subset_methods = base_methods
    else:
        all_methods = base_methods + [args.vinardo_method_name]
        common_subset_methods = all_methods

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
    methods_in_df = [
        col.split("_")[1]
        for col in results_df_top.columns
        if col.startswith("lddt_pli_") and col.split("_")[1] != "max" and col.split("_")[1] != "average"
    ]

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
