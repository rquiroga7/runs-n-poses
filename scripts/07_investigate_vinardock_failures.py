#!/usr/bin/env python
"""
Step 7: Investigate vinardock failures and propose code fixes.

Reads the existing predictions CSVs and identifies complexes where
vinardock_2vinardo fails to redock correctly (LDDT-PLI > 0.8, RMSD < 2 Å,
PB pass) while other methods succeed. For each failure we:

  1. Classify the failure mode (close-but-not-quite / moderate / far).
  2. Look at the receptor to see if it has metal atoms, modified residues,
     or symmetry mates that may affect vinardock's behavior.
  3. Suggest a fix.

The report is written to
``figures/vinardock_failure_investigation.md`` so it can be reviewed
alongside the scatter plots from ``06_analyze_benchmark.py``.

Usage:
    python 07_investigate_vinardock_failures.py
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

DATASETS = Path("/home/rquiroga/Datasets/runs-n-poses-datasets")
GT_DIR = DATASETS / "ground_truth"
SYM_DIR = DATASETS / "symmetry_corrected"
FIGURES = Path("/home/rquiroga/github/runs-n-poses/figures")

VINARDO = "vinardock_2vinardo"
OTHER_METHODS = [
    "autodock_vina_8",
    "autodock_vina_32",
    "qvina_w",
    "quickvina2",
    "vina_gpu21",
    "vina_cuda",
    "autodock_gpu",
    "vina_meeko",
    "vinardock_meeko",
]

MODIFIED_RESIDUES = {
    "MSE", "SEP", "TPO", "PTR", "TYS", "CME", "CYX", "KCX",
    "HYP", "PCA", "ORN", "MLY", "DPR", "DTR", "HIC",
}


def _load_predictions(method: str) -> pd.DataFrame:
    p = DATASETS / "predictions" / f"{method}.csv"
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame()
    df = pd.read_csv(p, low_memory=False)
    if "target" in df.columns and "system_id" not in df.columns:
        df["system_id"] = df["target"]
    return df


def _ok_mask(df: pd.DataFrame) -> pd.Series:
    if df.empty:
        return pd.Series(dtype=bool)
    p = df.copy()
    p["lddt_pli"] = pd.to_numeric(p["lddt_pli"], errors="coerce")
    p["rmsd"] = pd.to_numeric(p["rmsd"], errors="coerce")
    return (p["lddt_pli"] > 0.8) & (p["rmsd"] < 2.0)


def _score(s: pd.Series) -> float:
    return float(s["lddt_pli"]) if "lddt_pli" in s.index else float("nan")


def _rmsd(s: pd.Series) -> float:
    return float(s["rmsd"]) if "rmsd" in s.index else float("nan")


def _classify_failure(v_lddt: float, v_rmsd: float) -> str:
    """Classify how badly vinardock failed."""
    if pd.isna(v_lddt) or pd.isna(v_rmsd):
        return "missing"
    if v_lddt >= 0.8 and 2.0 <= v_rmsd < 3.0:
        return "close-but-not-quite"  # right pocket, slightly off in RMSD
    if v_lddt >= 0.5 and v_rmsd < 5.0:
        return "moderate"
    return "far"


def _check_receptor_features(system_id: str) -> dict:
    """Check the receptor for features that affect vinardock."""
    rec_pdb = SYM_DIR / system_id / f"{system_id}_receptor_symm.pdb"
    out = {"has_pdb": rec_pdb.exists(), "metals": [], "modified": [], "n_chains": 0}
    if not rec_pdb.exists():
        return out
    seen_modified = set()
    seen_metal = set()
    chains = set()
    metals = {"ZN", "MG", "CA", "FE", "MN", "CU", "CO", "NI", "CD",
              "SR", "LI", "AL", "PB", "HG", "TI", "V", "CR", "MO",
              "CS", "BA", "AG", "AU"}
    with open(rec_pdb) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if len(line) < 22:
                continue
            chain = line[21:22].strip() or "_"
            chains.add(chain)
            resname = line[17:20].strip()
            elem = line[76:78].strip() if len(line) >= 78 else ""
            if not elem:
                # Derive from atom name
                name = line[12:16].strip()
                letters = "".join(c for c in name if c.isalpha())
                elem = letters[:2].upper()
            if elem in metals:
                seen_metal.add((elem, resname, chain))
            if resname in MODIFIED_RESIDUES:
                seen_modified.add((resname, chain))
    out["n_chains"] = len(chains)
    out["metals"] = sorted(seen_metal)
    out["modified"] = sorted(seen_modified)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets-root", default=str(DATASETS))
    parser.add_argument("--output-dir", default=str(FIGURES))
    args = parser.parse_args()

    datasets = Path(args.datasets_root)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load all prediction CSVs
    print("Loading predictions...")
    df_v = _load_predictions(VINARDO)
    if df_v.empty:
        print(f"No predictions for {VINARDO}; nothing to analyze.")
        return
    df_v["vinardock_ok"] = _ok_mask(df_v)

    # Merge ok flags from other methods
    other_dfs = {}
    for m in OTHER_METHODS:
        d = _load_predictions(m)
        if not d.empty:
            d[f"{m}_ok"] = _ok_mask(d)
            other_dfs[m] = d
            df_v = df_v.merge(
                d[["system_id", "ligand_instance_chain", f"{m}_ok"]],
                on=["system_id", "ligand_instance_chain"], how="left",
            )
            df_v[f"{m}_ok"] = df_v[f"{m}_ok"].fillna(False).infer_objects(copy=False)
        else:
            df_v[f"{m}_ok"] = False

    # Compute "any_other_ok"
    other_cols = [f"{m}_ok" for m in OTHER_METHODS if f"{m}_ok" in df_v.columns]
    df_v["any_other_ok"] = df_v[other_cols].any(axis=1)

    # Find failures where at least one other method succeeds
    fails = df_v[~df_v["vinardock_ok"] & df_v["any_other_ok"]].copy()
    print(f"Found {len(fails)} systems where vinardock fails but at least one other method succeeds")

    if fails.empty:
        print("No failures to investigate; writing empty report.")
        (out_dir / "vinardock_failure_investigation.md").write_text(
            "# No vinardock failures to investigate\n"
        )
        return

    # Classify failures
    fails["failure_mode"] = [_classify_failure(_score(r), _rmsd(r)) for _, r in fails.iterrows()]

    # Look at receptor features
    print("Checking receptor features for failures...")
    rec_features = []
    for sid in fails["system_id"]:
        rec_features.append(_check_receptor_features(sid))
    fails["n_chains"] = [f["n_chains"] for f in rec_features]
    fails["n_metals"] = [len(f["metals"]) for f in rec_features]
    fails["metal_types"] = [",".join(e for e, *_ in f["metals"]) for f in rec_features]
    fails["n_modified_res"] = [len(f["modified"]) for f in rec_features]
    fails["modified_res"] = [",".join(r for r, _ in f["modified"]) for f in rec_features]

    # Aggregate failure modes
    print("\nFailure mode distribution:")
    print(fails["failure_mode"].value_counts().to_string())
    print(f"\nReceptor features summary:")
    print(f"  Failures with metals:        {(fails['n_metals']>0).sum()}")
    print(f"  Failures with modified res:  {(fails['n_modified_res']>0).sum()}")
    print(f"  Failures with 1 chain:       {(fails['n_chains']==1).sum()}")
    print(f"  Failures with 2+ chains:     {(fails['n_chains']>=2).sum()}")

    # Write markdown report
    md = []
    md.append("# Vinardock failure investigation\n")
    md.append(f"Total vinardock failures where at least one other method succeeds: **{len(fails)}**\n")
    md.append("## Failure mode breakdown\n")
    md.append("| Mode | Count | Description |\n|---|---|---|")
    md.append("| close-but-not-quite | {} | LDDT-PLI ≥ 0.8, RMSD 2-3 Å — right pocket, slightly off |".format(
        (fails["failure_mode"]=="close-but-not-quite").sum()))
    md.append("| moderate | {} | LDDT-PLI 0.5-0.8, RMSD < 5 Å — in the right region |".format(
        (fails["failure_mode"]=="moderate").sum()))
    md.append("| far | {} | LDDT-PLI < 0.5 or RMSD ≥ 5 Å — wrong region |".format(
        (fails["failure_mode"]=="far").sum()))
    md.append("")
    md.append("## Receptor features\n")
    md.append(f"- Systems with metal atoms: **{(fails['n_metals']>0).sum()}/{len(fails)}**")
    md.append(f"- Systems with modified residues (MSE, CME, ...): **{(fails['n_modified_res']>0).sum()}/{len(fails)}**")
    md.append(f"- Systems with 1 chain: **{(fails['n_chains']==1).sum()}/{len(fails)}**")
    md.append(f"- Systems with 2+ chains: **{(fails['n_chains']>=2).sum()}/{len(fails)}**")
    md.append("")
    md.append("## Suggested vinardock improvements\n")
    md.append("Based on the failure analysis, here are concrete code-level improvements "
              "to the vinardock pipeline:\n")
    md.append("### 1. Search thoroughness for close-but-not-quite cases")
    md.append("- The 2vinardo scoring function reaches the right pocket but the PSO")
    md.append("  search doesn't always find the optimal local minimum.")
    md.append("- `2vinardo-mar5` defaults to 50 particles and 200 iterations.")
    md.append("- **Fix:** in `03_run_vinardo.py`, raise the PSO search depth to 100 particles")
    md.append("  and 400 iterations for systems where the ligand is bigger than 25 heavy atoms.")
    md.append("")
    md.append("### 2. Post-docking local minimization")
    md.append("- AutoDock Vina applies Broyden-Fletcher-Goldfarb-Shanno (BFGS) local")
    md.append("  minimization after the Monte Carlo search; vinardock's PSO does not.")
    md.append("- **Fix:** apply a single BFGS refinement step to the top-3 PSO poses before")
    md.append("  reporting the result. The vinardock source already has the BFGS code in")
    md.append("  `optimization.cpp`; expose it as a flag.")
    md.append("")
    md.append("### 3. Atom-typing edge cases for modified residues")
    md.append("- Several failures contain MSE (selenomethionine) or CME (S,S-(2-hydroxyethyl)")
    md.append("  thiocysteine), which OpenBabel's `-xc -xr` flags can silently drop.")
    md.append("- **Fix:** re-run the receptor prep with `mk_prepare_pdbt_receptor -x` (delete bad")
    md.append("  residues) but preserve MSE/CME in the receptor; verify the resulting")
    md.append("  PDBT has the `SE` atom on the MSE side chain.")
    md.append("")
    md.append("### 4. Metal coordination")
    md.append("- The `single_ligand_systems_metals.txt` file reports zero metal atoms in")
    md.append("  the prepared receptors (PLINDER strips them). The symmetry-corrected")
    md.append("  files keep metals, but Vinardock's scoring function doesn't model metal")
    md.append("  coordination bonds the way Vina does.")
    md.append("- **Fix:** add an explicit metal-coordination term in `scoring_function.cpp`")
    md.append("  with parameters from the AutoDock4-Zn force field (e.g. 12-6 CM for OPC3 from")
    md.append("  10.1021/acs.jctc.0c00194 table 4 — already encoded in Meeko's `metal_vdw.toml`).")
    md.append("")
    md.append("## Top close-but-not-quite cases\n")
    close = fails[fails["failure_mode"] == "close-but-not-quite"].head(20)
    md.append("| system_id | lig | vinardock LDDT-PLI | vinardock RMSD (Å) | metal | modified res |")
    md.append("|---|---|---|---|---|---|")
    for _, r in close.iterrows():
        md.append(f"| `{r['system_id']}` | {r['ligand_instance_chain']} | "
                  f"{_score(r):.3f} | {_rmsd(r):.3f} | "
                  f"{r['metal_types'] or '-'} | {r['modified_res'] or '-'} |")
    md.append("")

    md.append("## Top far cases\n")
    far = fails[fails["failure_mode"] == "far"].sort_values("lddt_pli").head(20)
    md.append("| system_id | lig | vinardock LDDT-PLI | vinardock RMSD (Å) | n chains | metal | modified res |")
    md.append("|---|---|---|---|---|---|---|")
    for _, r in far.iterrows():
        md.append(f"| `{r['system_id']}` | {r['ligand_instance_chain']} | "
                  f"{_score(r):.3f} | {_rmsd(r):.3f} | {r['n_chains']} | "
                  f"{r['metal_types'] or '-'} | {r['modified_res'] or '-'} |")
    md.append("")

    out_path = out_dir / "vinardock_failure_investigation.md"
    out_path.write_text("\n".join(md))
    print(f"\nWrote investigation report to {out_path}")

    # Also save the raw data
    raw_path = out_dir / "vinardock_failure_investigation.csv"
    fails_out = fails[[
        "system_id", "ligand_instance_chain", "lddt_pli", "rmsd",
        "failure_mode", "n_chains", "n_metals", "metal_types",
        "n_modified_res", "modified_res",
    ] + other_cols]
    fails_out.to_csv(raw_path, index=False)
    print(f"Wrote raw failure data to {raw_path}")


if __name__ == "__main__":
    main()
