#!/usr/bin/env python
"""
Step 6: Aggregate the docking benchmark into per-method speed / accuracy
metrics and produce a ``speedup vs redocking%`` scatter plot.

For each method listed in ``--methods`` the script reads:

  * ``predictions/<method>.csv`` for accuracy metrics
    (``lddt_pli``, ``rmsd``).
  * ``posebusters_results/<method>.csv`` for ``pb_success``.
  * ``runtime.json`` written into every ``<output>/<sys>_<lig>/`` folder
    by the docking runners.

The redocking success rate is the fraction of systems that satisfy all
of:

  * ``lddt_pli > 0.8``
  * ``rmsd < 2.0``
  * PoseBusters ``pb_success == 1``

The "speedup" is computed relative to ``autodock_vina_8`` (the default
reference) using the mean per-system runtime.

Outputs:

  * ``benchmark_summary.csv``: per-method success / runtime table
  * ``benchmark_speedup_vs_redocking.png``: log-scale speedup vs %
    redocking scatter
  * ``vinardock_failure_analysis.csv``: per-system reasons why vinardock
    fails where other methods succeed (used for the "possible fixes"
    investigation step).

Usage:
    python 06_analyze_benchmark.py \\
        --datasets-root /home/rquiroga/Datasets/runs-n-poses-datasets \\
        --output-dir /home/rquiroga/github/runs-n-poses/figures \\
        --reference autodock_vina_8
"""

import argparse
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DATASETS_DEFAULT = Path("/home/rquiroga/Datasets/runs-n-poses-datasets")
FIGURES_DEFAULT = Path("/home/rquiroga/github/runs-n-poses/figures")

DEFAULT_METHODS = [
    "autodock_vina_8",
    "autodock_vina_32",
    "qvina_w",
    "quickvina2",
    "autodock_gpu",
    "autodock_gpu_mostpop",
    "vina_8_meeko",
    "vina_32_meeko",
    "vinardock_2vinardo",
    "vinardock_meeko",
    "rdock",
    "gnina",
]

METHOD_DISPLAY = {
    "autodock_vina_8": "AutoDock Vina (exh=8)",
    "autodock_vina_32": "AutoDock Vina (exh=32)",
    "qvina_w": "QuickVina-W",
    "quickvina2": "QuickVina2",
    "autodock_gpu": "AutoDock-GPU (best energy)",
    "autodock_gpu_mostpop": "AutoDock-GPU (most pop. cluster)",
    "vina_8_meeko": "Vina (Meeko prep, exh=8)",
    "vina_32_meeko": "Vina (Meeko prep, exh=32)",
    "vinardock_2vinardo": "Vinardo (obabel prep)",
    "vinardock_meeko": "Vinardo (Meeko prep)",
    "rdock": "rDock",
    "gnina": "GNINA",
}

METHOD_COLORS = {
    "autodock_vina_8": "#1f77b4",
    "autodock_vina_32": "#aec7e8",
    "qvina_w": "#ff7f0e",
    "quickvina2": "#ffbb78",
    "autodock_gpu": "#d62728",
    "autodock_gpu_mostpop": "#ff9896",
    "vina_8_meeko": "#9467bd",
    "vina_32_meeko": "#c5b0d5",
    "vinardock_2vinardo": "#8c564b",
    "vinardock_meeko": "#e377c2",
    "rdock": "#17becf",
    "gnina": "#bcbd22",
}

METHOD_MARKERS = {
    "autodock_vina_8": "o",
    "autodock_vina_32": "o",
    "qvina_w": "s",
    "quickvina2": "D",
    "autodock_gpu": "P",
    "autodock_gpu_mostpop": "p",
    "vina_8_meeko": "X",
    "vina_32_meeko": "d",
    "vinardock_2vinardo": "*",
    "vinardock_meeko": "h",
    "rdock": ">",
    "gnina": "1",
}


METHOD_ALIASES = {
    "autodock_vina_8": "vina",
}

def _load_predictions(datasets: Path, method: str) -> pd.DataFrame:
    csv_name = "autodock_gpu.csv" if method == "autodock_gpu_mostpop" else f"{method}.csv"
    pred_csv = datasets / "predictions" / csv_name
    if not pred_csv.exists() or pred_csv.stat().st_size == 0:
        return pd.DataFrame()
    try:
        df = pd.read_csv(pred_csv, low_memory=False)
        if "target" in df.columns and "system_id" not in df.columns:
            df["system_id"] = df["target"]
        # Filter to this method's rows if CSV mixes multiple methods
        if "method" in df.columns:
            methods_in_csv = set(df["method"].unique())
            if method in methods_in_csv:
                df = df[df["method"] == method]
            # Also include known alias rows (e.g. 'vina' for 'autodock_vina_8')
            alias = METHOD_ALIASES.get(method)
            if alias and alias in methods_in_csv:
                alias_rows = pd.read_csv(pred_csv, low_memory=False)
                if "target" in alias_rows.columns and "system_id" not in alias_rows.columns:
                    alias_rows["system_id"] = alias_rows["target"]
                alias_rows = alias_rows[alias_rows["method"] == alias]
                df = pd.concat([df, alias_rows], ignore_index=True)
        return df
    except Exception:
        return pd.DataFrame()


def _load_posebusters(datasets: Path, method: str) -> pd.DataFrame:
    pb_csv = datasets / "posebusters_results" / f"{method}.csv"
    if not pb_csv.exists() or pb_csv.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(pb_csv, low_memory=False)
    except Exception:
        return pd.DataFrame()


def _gather_runtimes(datasets: Path, method: str) -> pd.DataFrame:
    """Walk the per-method output dir and collect runtime.json files."""
    method_dir_candidates = {
        "vinardock_2vinardo": datasets / "vinardo_outputs",
        "vinardock_meeko": datasets / "vinardock_meeko",
        "vina_meeko": datasets / "vina_meeko",
        "autodock_gpu_mostpop": datasets / "autodock_gpu",
        "rdock": datasets / "rdock",
        "gnina": datasets / "gnina",
    }
    method_dir = method_dir_candidates.get(method, datasets / method)
    if not method_dir.is_dir():
        return pd.DataFrame()

    search_root = method_dir

    rows = []
    if not search_root.is_dir():
        return pd.DataFrame(rows)

    # Check flat layout: <sys_id>/runtime.json (autodock_gpu, vina_meeko, etc.)
    flat_count = 0
    for d in search_root.iterdir():
        if not d.is_dir():
            continue
        rt_file = d / "runtime.json"
        if rt_file.exists():
            flat_count += 1
            try:
                with open(rt_file) as f:
                    rows.append({"system_id": d.name,
                                 "runtime_seconds": json.load(f).get("runtime_seconds")})
            except Exception:
                pass

    # If few flat entries, also check nested: <sys_id>/<sys_id>_<lig>/runtime.json (vina-style)
    if flat_count < len(list(search_root.iterdir())) // 2:
        for d in search_root.iterdir():
            if not d.is_dir():
                continue
            for sub in d.iterdir():
                if not sub.is_dir():
                    continue
                rt_file = sub / "runtime.json"
                if rt_file.exists():
                    # Deduplicate: skip if sys_id already found via flat layout
                    if any(r["system_id"] == d.name for r in rows):
                        continue
                    try:
                        with open(rt_file) as f:
                            rows.append({"system_id": d.name,
                                         "ligand_chain": sub.name.split("_")[-1],
                                         "runtime_seconds": json.load(f).get("runtime_seconds")})
                    except Exception:
                        pass
    return pd.DataFrame(rows)


def _success_mask(df_pred: pd.DataFrame, df_pb: pd.DataFrame) -> pd.Series:
    """Boolean Series: True iff system is correctly redocked."""
    if df_pred.empty:
        return pd.Series(dtype=bool)
    p = df_pred.copy()
    p["lddt_pli"] = pd.to_numeric(p["lddt_pli"], errors="coerce")
    p["rmsd"] = pd.to_numeric(p["rmsd"], errors="coerce")
    if not df_pb.empty and "pb_success" in df_pb.columns:
        pb = df_pb.copy()
        pb["pb_success"] = pd.to_numeric(pb["pb_success"], errors="coerce")
        if "system_id" in pb.columns and "ligand_instance_chain" in pb.columns:
            pb = pb.drop_duplicates(subset=["system_id", "ligand_instance_chain"])
            p = p.merge(
                pb[["system_id", "ligand_instance_chain", "pb_success"]],
                on=["system_id", "ligand_instance_chain"],
                how="left",
            )
    mask = (p["lddt_pli"] > 0.8) & (p["rmsd"] < 2.0)
    if "pb_success" in p.columns:
        pb_mask = p["pb_success"].fillna(1).astype(float)
        mask = mask & (pb_mask == 1.0)
    return mask


def compute_summary(datasets: Path, methods: list, annotations: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for m in methods:
        df_pred = _load_predictions(datasets, m)
        df_pb = _load_posebusters(datasets, m)
        df_rt = _gather_runtimes(datasets, m)

        n_pred = df_pred["system_id"].nunique() if not df_pred.empty else 0
        mask = _success_mask(df_pred, df_pb)
        # Count unique systems where at least one row passes the mask (top-1)
        n_success = df_pred.loc[mask, "system_id"].nunique() if not df_pred.empty and not mask.empty else 0
        success_pct = 100.0 * n_success / max(1, n_pred) if n_pred else float("nan")
        mean_rt = df_rt["runtime_seconds"].dropna().mean() if not df_rt.empty else float("nan")
        median_rt = df_rt["runtime_seconds"].dropna().median() if not df_rt.empty else float("nan")

        rows.append({
            "method": m,
            "display": METHOD_DISPLAY.get(m, m),
            "n_predicted": n_pred,
            "n_success": n_success,
            "success_pct": success_pct,
            "mean_runtime_s": mean_rt,
            "median_runtime_s": median_rt,
        })
    return pd.DataFrame(rows)


def plot_speedup_vs_redocking(summary: pd.DataFrame, reference: str, output_path: Path):
    if summary.empty:
        print("Summary is empty; nothing to plot.")
        return

    ref_row = summary[summary["method"] == reference]
    if ref_row.empty:
        print(f"Reference method {reference} not in summary; cannot compute speedup.")
        return
    ref_rt = ref_row["mean_runtime_s"].iloc[0]
    if not np.isfinite(ref_rt) or ref_rt == 0:
        print("Reference runtime is invalid; cannot compute speedup.")
        return

    fig, ax = plt.subplots(figsize=(8, 6))
    for _, row in summary.iterrows():
        if not np.isfinite(row["mean_runtime_s"]) or not np.isfinite(row["success_pct"]):
            continue
        speedup = ref_rt / row["mean_runtime_s"]
        ax.scatter(
            speedup, row["success_pct"],
            color=METHOD_COLORS.get(row["method"], "#333333"),
            marker=METHOD_MARKERS.get(row["method"], "o"),
            s=140, edgecolor="black", linewidth=0.7,
            label=METHOD_DISPLAY.get(row["method"], row["method"]),
        )
        ax.annotate(
            METHOD_DISPLAY.get(row["method"], row["method"]),
            (speedup, row["success_pct"]),
            textcoords="offset points", xytext=(8, 5),
            fontsize=8, alpha=0.9,
        )
    ax.set_xscale("log")
    ax.set_xlabel(f"Speedup vs {METHOD_DISPLAY.get(reference, reference)} (log)")
    ax.set_ylabel("Redocking success (%)  [LDDT-PLI>0.8, RMSD<2 Å, PB pass]")
    ax.set_title("Runs N' Poses docking benchmark: speedup vs accuracy")
    ax.grid(True, which="both", ls="--", alpha=0.3)
    ax.axhline(100, color="gray", ls=":", alpha=0.5)
    ax.axvline(1, color="gray", ls=":", alpha=0.5)
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Saved speedup scatter to {output_path}")


def vinardock_failure_analysis(datasets: Path, output_path: Path,
                                vinardo_method: str = "vinardock_2vinardo",
                                other_methods: list | None = None) -> pd.DataFrame:
    """Identify complexes where vinardock fails but other methods succeed.

    Output: per-system, per-other-method breakdown of which checks vinardock
    missed compared to the rest of the field. Useful for the
    ``look for possible fixes`` step.
    """
    other_methods = other_methods or [
        "autodock_vina_8", "autodock_vina_32", "qvina_w", "quickvina2",
        "vina_gpu21", "vina_cuda", "autodock_gpu", "rdock", "gnina",
    ]
    pred_v = _load_predictions(datasets, vinardo_method)
    if pred_v.empty:
        print(f"No predictions for {vinardo_method}; skipping failure analysis.")
        return pd.DataFrame()

    pred_v = pred_v.copy()
    pred_v["lddt_pli"] = pd.to_numeric(pred_v["lddt_pli"], errors="coerce")
    pred_v["rmsd"] = pd.to_numeric(pred_v["rmsd"], errors="coerce")
    pred_v["vinardo_ok"] = (pred_v["lddt_pli"] > 0.8) & (pred_v["rmsd"] < 2.0)

    rows = []
    for om in other_methods:
        df = _load_predictions(datasets, om)
        if df.empty:
            continue
        df = df.copy()
        df["lddt_pli"] = pd.to_numeric(df["lddt_pli"], errors="coerce")
        df["rmsd"] = pd.to_numeric(df["rmsd"], errors="coerce")
        df[f"{om}_ok"] = (df["lddt_pli"] > 0.8) & (df["rmsd"] < 2.0)
        keep = ["system_id", "ligand_instance_chain", f"{om}_ok"]
        pred_v = pred_v.merge(df[keep],
                              on=["system_id", "ligand_instance_chain"],
                              how="left")
        pred_v[f"{om}_ok"] = pred_v[f"{om}_ok"].fillna(False).infer_objects(copy=False)

    any_other_ok = pred_v[[f"{m}_ok" for m in other_methods if f"{m}_ok" in pred_v.columns]].any(axis=1)
    fails = pred_v[~pred_v["vinardo_ok"] & any_other_ok].copy()
    print(f"Found {len(fails)} systems where {vinardo_method} fails but another method succeeds.")

    if fails.empty:
        return fails

    # Compute per-system aggregate
    fail_table = fails[["system_id", "ligand_instance_chain", "lddt_pli", "rmsd"]].copy()
    for om in other_methods:
        col = f"{om}_ok"
        if col in fails.columns:
            fail_table[f"{om}_ok"] = fails[col].astype(int)
    fail_table["n_other_ok"] = fail_table[[f"{m}_ok" for m in other_methods
                                            if f"{m}_ok" in fail_table.columns]].sum(axis=1)

    # Categorize the failure mode using vina_8 RMSD and LDDT-PLI deltas
    if not _load_predictions(datasets, "autodock_vina_8").empty:
        ref = _load_predictions(datasets, "autodock_vina_8").copy()
        ref["lddt_pli"] = pd.to_numeric(ref["lddt_pli"], errors="coerce")
        ref["rmsd"] = pd.to_numeric(ref["rmsd"], errors="coerce")
        ref = ref.rename(columns={"lddt_pli": "ref_lddt_pli", "rmsd": "ref_rmsd"})
        fail_table = fail_table.merge(
            ref[["system_id", "ligand_instance_chain", "ref_lddt_pli", "ref_rmsd"]],
            on=["system_id", "ligand_instance_chain"], how="left",
        )
        fail_table["lddt_pli_gap"] = fail_table["ref_lddt_pli"] - fail_table["lddt_pli"]
        fail_table["rmsd_gap"] = fail_table["rmsd"] - fail_table["ref_rmsd"]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fail_table.to_csv(output_path, index=False)
    print(f"Saved failure analysis to {output_path}")
    return fail_table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets-root", type=str, default=str(DATASETS_DEFAULT))
    parser.add_argument("--output-dir", type=str, default=str(FIGURES_DEFAULT))
    parser.add_argument("--methods", type=str, default=",".join(DEFAULT_METHODS),
                        help="comma-separated list of method names to include")
    parser.add_argument("--reference", type=str, default="autodock_vina_8",
                        help="reference method for speedup calculation")
    parser.add_argument("--system-list", type=str, default=None,
                        help="optional system-list to restrict the analysis")
    args = parser.parse_args()

    datasets = Path(args.datasets_root)
    out_dir = Path(args.output_dir)
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]

    annotations = pd.read_csv(datasets / "annotations.csv")
    annotations["release_date"] = pd.to_datetime(annotations["release_date"])
    if args.system_list:
        with open(args.system_list) as f:
            allowed = set(l.strip() for l in f if l.strip())
        annotations = annotations[annotations["system_id"].astype(str).isin(allowed)].reset_index(drop=True)

    print("Computing per-method summary...")
    summary = compute_summary(datasets, methods, annotations)
    summary_path = out_dir / "benchmark_summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"Saved summary to {summary_path}")
    print(summary.to_string(index=False))

    # Speedup vs redocking scatter
    plot_speedup_vs_redocking(summary, args.reference,
                              out_dir / "benchmark_speedup_vs_redocking.png")

    # Vinardock failure analysis
    vinardock_failure_analysis(datasets, out_dir / "vinardock_failure_analysis.csv")

    # Also produce a bar chart of redocking %
    fig, ax = plt.subplots(figsize=(10, 5))
    order = summary.sort_values("success_pct", ascending=False)
    ax.barh([METHOD_DISPLAY.get(m, m) for m in order["method"]], order["success_pct"],
            color=[METHOD_COLORS.get(m, "#333333") for m in order["method"]])
    ax.set_xlabel("Redocking success (%)  [LDDT-PLI>0.8, RMSD<2 Å, PB pass]")
    ax.set_title("Per-method redocking success on Runs N' Poses (single ligand)")
    ax.invert_yaxis()
    plt.tight_layout()
    out_path = out_dir / "benchmark_redocking_bars.png"
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"Saved redocking bar chart to {out_path}")

    # Runtime log-scale bar
    fig, ax = plt.subplots(figsize=(10, 5))
    rt = summary.copy()
    rt = rt.dropna(subset=["mean_runtime_s"]).sort_values("mean_runtime_s")
    ax.barh([METHOD_DISPLAY.get(m, m) for m in rt["method"]], rt["mean_runtime_s"],
            color=[METHOD_COLORS.get(m, "#333333") for m in rt["method"]])
    ax.set_xscale("log")
    ax.set_xlabel("Mean per-system runtime (s, log scale)")
    ax.set_title("Docking method runtimes (mean per system, log scale)")
    plt.tight_layout()
    out_path = out_dir / "benchmark_runtime_bars.png"
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"Saved runtime bar chart to {out_path}")


if __name__ == "__main__":
    main()
