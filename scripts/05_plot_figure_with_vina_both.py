#!/usr/bin/env python
"""
Step 5: Produce Figure 1 including both `autodock_vina_8` and `autodock_vina_32`.

This script loads predictions and PoseBusters outputs for Vinardock, both Vina runs,
and AF3, and calls the shared plotting utilities.
"""

import os
os.environ["MPLBACKEND"] = "Agg"

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent.parent))
import plotting


def load_data(data_dir: Path, methods: list):
    annotated_df = pd.read_csv(data_dir / "annotations.csv")
    annotated_df["release_date"] = pd.to_datetime(annotated_df["release_date"])

    bust_dfs = {}
    for m in methods:
        filename = data_dir / "posebusters_results" / f"{m}.csv"
        if filename.exists():
            # Some PoseBusters CSVs may contain malformed rows (extra commas).
            # Use the python engine and skip bad lines to avoid crashing the plot loader.
            try:
                bust_df = pd.read_csv(filename, low_memory=False, dtype=str, on_bad_lines='skip')
            except TypeError:
                # Older pandas versions use error_bad_lines; try fallback
                bust_df = pd.read_csv(filename, low_memory=False, dtype=str, error_bad_lines=False)
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
            "runtime_seconds",
            "exhaustiveness",
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
        # Ensure string types before concatenation to avoid dtype errors
        full_datasets[method]["group_key"] = (
            full_datasets[method]["system_id"].fillna("")
            + "__"
            + full_datasets[method]["ligand_instance_chain"].fillna("")
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
    # Pivot method-specific columns into wide format
    pivoted = df.pivot(
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
            "runtime_seconds",
            "exhaustiveness",
            "pred_pocket_f1",
        ],
    ).reset_index()
    pivoted.columns = [f"{col[0]}_{col[1]}" if len(col[1]) else col[0] for col in pivoted.columns]

    # Ensure annotated_df contains `group_key` for merging
    ann = annotated_df.copy()
    if "group_key" not in ann.columns:
        ann["group_key"] = (
            ann["system_id"].astype(str).fillna("")
            + "__"
            + ann["ligand_instance_chain"].astype(str).fillna("")
        )

    # Merge pivoted results into the full annotation table so that all annotated
    # single-ligand systems are present even if a method has no prediction for them.
    merged = ann.merge(pivoted, on="group_key", how="left")

    if "ligand_is_proper" not in merged.columns:
        proper_cols = [c for c in ["ligand_is_proper_x", "ligand_is_proper_y"] if c in merged.columns]
        if proper_cols:
            merged["ligand_is_proper"] = merged[proper_cols[0]]
            for col in proper_cols[1:]:
                merged["ligand_is_proper"] = merged["ligand_is_proper"].fillna(merged[col])
        else:
            merged["ligand_is_proper"] = False

    # Keep only proper ligands for plotting
    merged = merged[merged["ligand_is_proper"].fillna(False)].reset_index(drop=True)
    return merged


def main():
    parser = argparse.ArgumentParser(description="Plot Figure 1 including both Vina runs")
    parser.add_argument(
        "--data-dir",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets",
    )
    parser.add_argument(
        "--system-list",
        "--system_list",
        dest="system_list",
        type=str,
        default=None,
        help="Optional path to a file with one system_id per line to restrict plotting",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="/home/rquiroga/github/runs-n-poses/figures/fig_1_with_vina_both.png",
    )
    parser.add_argument(
        "--vina-method-names",
        type=str,
        default="autodock_vina_8,autodock_vina_32",
        help="Comma-separated method names for the two Vina runs",
    )
    parser.add_argument(
        "--vinardo-method-name",
        type=str,
        default="vinardock_2vinardo",
    )
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    # optional system list to restrict plotting to a subset of systems
    allowed_systems = None
    if args.system_list:
        slf = Path(args.system_list)
        if slf.exists():
            allowed_systems = set([l.strip() for l in slf.read_text().splitlines() if l.strip()])
            print(f"Restricting plot to {len(allowed_systems)} systems from: {slf}")
        else:
            print(f"Warning: --system-list file not found: {slf}; ignoring")
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    vina_names = [n.strip() for n in args.vina_method_names.split(",") if n.strip()]
    vinardo_name = args.vinardo_method_name
    af3_name = "af3"

    base_methods = [vinardo_name] + vina_names + [af3_name]
    all_methods = base_methods
    common_subset_methods = base_methods

    print(f"Loading data for methods: {all_methods}")

    full_datasets, annotated_df = load_data(data_dir, all_methods)

    # Restrict to single-ligand systems only (plot expects single-ligand complexes)
    if "num_ligand_chains" in annotated_df.columns:
        before = len(annotated_df)
        annotated_df = annotated_df[annotated_df["num_ligand_chains"] == 1].reset_index(drop=True)
        print(f"Filtered annotations to single-ligand systems: {before} -> {len(annotated_df)}")
    # apply system-list filtering to annotations and loaded prediction/pb data
    if allowed_systems is not None:
        annotated_df = annotated_df[annotated_df["system_id"].astype(str).isin(allowed_systems)].reset_index(drop=True)
        for m in list(full_datasets.keys()):
            df = full_datasets[m]
            if "system_id" in df.columns:
                df = df[df["system_id"].astype(str).isin(allowed_systems)].reset_index(drop=True)
            full_datasets[m] = df

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

    plotting.COLORS[vinardo_name] = "#2ca02c"
    plotting.SHAPES[vinardo_name] = "^"
    plotting.NAME_MAPPING[vinardo_name] = "Vinardock"

    # assign distinct colors for the two Vina runs
    if len(vina_names) >= 1:
        plotting.COLORS[vina_names[0]] = "#ff7f0e"
        plotting.SHAPES[vina_names[0]] = "v"
        suffix0 = vina_names[0].split('_')[-1]
        plotting.NAME_MAPPING[vina_names[0]] = f"Autodock Vina (exhaustiveness_{suffix0})"
    if len(vina_names) >= 2:
        plotting.COLORS[vina_names[1]] = "#d62728"
        plotting.SHAPES[vina_names[1]] = "s"
        suffix1 = vina_names[1].split('_')[-1]
        plotting.NAME_MAPPING[vina_names[1]] = f"Autodock Vina (exhaustiveness_{suffix1})"

    plotting.COLORS[af3_name] = plotting.COLORS.get(af3_name, plotting.COLORS.get('af3', '#1f77b4'))
    plotting.SHAPES[af3_name] = plotting.SHAPES.get(af3_name, 'o')

    print("\nCreating figure using plotting.make_main_figure()...")

    plotting.make_main_figure(
        results_df_top,
        str(output_path),
        methods=methods_in_df,
    )
    # --- New: compute per-method success % and average runtime; create scatter + table
    print("\nComputing runtime and success summary for methods...")
    summary = []
    # success condition: LDDT-PLI > 0.8, RMSD < 2.0, PB-valid == 1
    for m in methods_in_df:
        lddt_col = f"lddt_pli_{m}"
        rmsd_col = f"rmsd_{m}"
        pb_col = f"pb_success_{m}"
        rt_col = f"runtime_seconds_{m}"

        if lddt_col in results_df_top.columns and rmsd_col in results_df_top.columns and pb_col in results_df_top.columns:
            cond = (
                (results_df_top[lddt_col].astype(float) > 0.8)
                & (results_df_top[rmsd_col].astype(float) < 2.0)
                & (results_df_top[pb_col].astype(float) == 1.0)
            )
            success_pct = 100.0 * cond.sum() / max(1, len(results_df_top))
        else:
            success_pct = float('nan')

        # average runtime (seconds) — prefer pivoted column, fall back to full_datasets mean
        avg_rt = None
        if rt_col in results_df_top.columns:
            try:
                avg_rt = pd.to_numeric(results_df_top[rt_col], errors='coerce').dropna().mean()
            except Exception:
                avg_rt = None
        if avg_rt is None or np.isnan(avg_rt):
            # fallback: check full_datasets (if available)
            try:
                if m in full_datasets and "runtime_seconds" in full_datasets[m].columns:
                    avg_rt = pd.to_numeric(full_datasets[m]["runtime_seconds"], errors='coerce').dropna().mean()
            except Exception:
                avg_rt = None

        summary.append({"method": m, "success_pct": success_pct, "avg_runtime_s": float(avg_rt) if avg_rt is not None and not np.isnan(avg_rt) else pd.NA})

    summary_df = pd.DataFrame(summary)
    table_path = output_path.parent / "runtimes_success_table.csv"
    summary_df.to_csv(table_path, index=False)
    print(f"Saved runtime-success table to: {table_path}")

    # Scatter plot: x = success_pct, y = avg_runtime_s (log scale)
    scatter_path = output_path.parent / (output_path.stem + "_runtime_scatter.png")
    plt.figure(figsize=(6, 5))
    for _, row in summary_df.iterrows():
        m = row['method']
        x = row['success_pct']
        y = row['avg_runtime_s']
        if pd.isna(x) or pd.isna(y):
            continue
        label = plotting.NAME_MAPPING.get(m, m)
        color = plotting.COLORS.get(m, '#333333')
        marker = plotting.SHAPES.get(m, 'o')
        plt.scatter(x, y, label=label, color=color, marker=marker, s=80)
        plt.text(x, y, ' ' + label, verticalalignment='center', fontsize=9)
    plt.xlabel('Success % (LDDT-PLI>0.8, RMSD<2Å, PB pass)')
    plt.ylabel('Average runtime (s)')
    plt.yscale('log')
    plt.grid(True, which='both', ls='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(scatter_path, dpi=200)
    plt.close()
    print(f"Saved runtime-success scatter to: {scatter_path}")

    print(f"\nFigure saved to: {output_path}")
    print("Done!")


if __name__ == "__main__":
    main()
