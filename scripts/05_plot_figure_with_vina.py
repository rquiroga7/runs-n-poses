#!/usr/bin/env python
"""
Step 5: Reproduce Figure 1 replacing Boltz-1 with AutoDock Vina.

This mirrors `05_plot_figure_with_vinardo.py` but substitutes the `boltz` method
with `vina` (predictions/vina.csv and posebusters_results/vina.csv).
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
    annotated_df = pd.read_csv(data_dir / "annotations.csv")
    annotated_df["release_date"] = pd.to_datetime(annotated_df["release_date"])

    bust_dfs = {}
    for m in methods:
        filename = data_dir / "posebusters_results" / f"{m}.csv"
        if filename.exists():
            bust_df = pd.read_csv(filename)
            if "ligand_chain" in bust_df.columns and "ligand_instance_chain" not in bust_df.columns:
                bust_df["ligand_instance_chain"] = bust_df["ligand_chain"]
            checks = [
                'sanitization', 'inchi_convertible', 'all_atoms_connected',
                'bond_lengths', 'bond_angles', 'internal_steric_clash',
                'aromatic_ring_flatness', 'double_bond_flatness',
                'internal_energy', 'protein-ligand_maximum_distance',
                'minimum_distance_to_protein'
            ]
            if "pb_success" not in bust_df.columns:
                available = [c for c in checks if c in bust_df.columns]
                if available:
                    temp = bust_df[available].replace({'True': True, 'False': False})
                    bust_df["pb_success"] = temp.eq(True).all(axis=1).astype(float)
                else:
                    bust_df["pb_success"] = -1.0
            bust_dfs[m] = bust_df

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
        full_datasets[method] = df[keep_columns].rename(columns={"target": "system_id"}).reset_index(drop=True)
        full_datasets[method]["system_id"] = full_datasets[method]["system_id"].astype(str)
        full_datasets[method]["ligand_instance_chain"] = full_datasets[method]["ligand_instance_chain"].astype(str)
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
            bust_dfs[method]["system_id"] = bust_dfs[method]["system_id"].astype(str)
            bust_dfs[method]["ligand_instance_chain"] = bust_dfs[method]["ligand_instance_chain"].astype(str)
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
        description="Reproduce Figure 1 with AutoDock Vina replacing Boltz-1"
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
        default="/home/rquiroga/github/runs-n-poses/figures/fig_1_with_vina.png",
        help="Path to output figure file",
    )
    parser.add_argument(
        "--vina-method-name",
        type=str,
        default="vina",
        help="Method name for vina in predictions",
    )
    parser.add_argument(
        "--vinardo-method-name",
        type=str,
        default="vinardock_2vinardo",
        help="Method name for vinardock in predictions",
    )
    parser.add_argument(
        "--common-subset-only",
        action="store_true",
        help="Only include systems where ALL methods (including vina) have results",
    )

    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # We want to plot exactly: Vinardock, Autodock Vina, and AlphaFold3
    vinardo_name = args.vinardo_method_name
    vina_name = args.vina_method_name
    af3_name = "af3"

    # Order: Vinardock, Autodock Vina, AlphaFold3
    base_methods = [vinardo_name, vina_name, af3_name]
    all_methods = base_methods
    common_subset_methods = base_methods

    print(f"Loading data for methods: {all_methods}")
    print(f"Common subset methods: {common_subset_methods}")

    full_datasets, annotated_df = load_data(data_dir, all_methods)

    if not full_datasets:
        print("Error: No prediction data loaded")
        sys.exit(1)

    print(f"Loaded {len(full_datasets)} methods: {list(full_datasets.keys())}")

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

    prefix = "lddt_pli_"
    methods_in_df = list(set([
        col[len(prefix):]
        for col in results_df_top.columns
        if col.startswith(prefix) and col[len(prefix):] not in ("max", "average")
    ]))
    methods_in_df = [m for m in base_methods if m in methods_in_df]
    common_subset_methods = [m for m in base_methods if f"lddt_pli_{m}" in results_df_top.columns]

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

    # Assign colors/markers and display names for the three methods
    plotting.COLORS[vinardo_name] = "#2ca02c"
    plotting.SHAPES[vinardo_name] = "^"
    plotting.NAME_MAPPING[vinardo_name] = "Vinardock"

    plotting.COLORS[vina_name] = "#ff7f0e"
    plotting.SHAPES[vina_name] = "v"
    plotting.NAME_MAPPING[vina_name] = "Autodock Vina"

    # `af3` label already exists, but ensure shape/color are present
    plotting.COLORS[af3_name] = plotting.COLORS.get(af3_name, plotting.COLORS.get('af3', '#1f77b4'))
    plotting.SHAPES[af3_name] = plotting.SHAPES.get(af3_name, 'o')

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
