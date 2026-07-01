#!/usr/bin/env python
"""
Step 1: Generate the comprehensive system list for symmetry-corrected docking.

Produces two files from annotations.csv:
  - systems_for_symmetry_docking.txt   — one system_id per line
  - systems_for_symmetry_docking.csv   — system_id,proper_ligand_chain

Selection criteria:
  - Single-ligand systems (num_ligand_chains == 1)
  - Multi-ligand systems with exactly 1 proper (drug-like) ligand
"""

import csv
import os
from collections import defaultdict

ANNOTATIONS = os.path.join(
    os.path.dirname(__file__),
    "..",
    "runs-n-poses-datasets",
    "annotations.csv",
)
OUT_TXT = os.path.join(os.path.dirname(__file__), "systems_for_symmetry_docking.txt")
OUT_CSV = os.path.join(os.path.dirname(__file__), "systems_for_symmetry_docking.csv")


def main():
    with open(ANNOTATIONS) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    by_system = defaultdict(list)
    for r in rows:
        by_system[r["system_id"]].append(r)

    systems = []  # (system_id, proper_ligand_chain)
    for sys_id, ligands in sorted(by_system.items()):
        num_chains = int(ligands[0].get("num_ligand_chains", len(ligands)))
        proper = [l for l in ligands if l["ligand_is_proper"] == "True"]

        if num_chains == 1 and len(proper) == 1:
            # Single-ligand: the only chain is the proper one
            systems.append((sys_id, proper[0]["ligand_instance_chain"]))
        elif num_chains > 1 and len(proper) == 1:
            # Multi-ligand with exactly 1 proper ligand
            systems.append((sys_id, proper[0]["ligand_instance_chain"]))
        # Otherwise skip (no proper, or >1 proper)

    print(f"Total systems selected: {len(systems)}")

    with open(OUT_TXT, "w") as f:
        for sys_id, _ in systems:
            f.write(sys_id + "\n")

    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["system_id", "proper_ligand_chain"])
        for sys_id, chain in systems:
            writer.writerow([sys_id, chain])

    print(f"Wrote {OUT_TXT} ({len(systems)} lines)")
    print(f"Wrote {OUT_CSV} ({len(systems)} data rows)")


if __name__ == "__main__":
    main()
