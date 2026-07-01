#!/usr/bin/env python
"""Compute docking box from ground-truth ligand SDF: center + size (+10Å margin)."""
import math
from pathlib import Path


def parse_ligand_coords(lig_sdf: str | Path) -> list[tuple[float, float, float]]:
    """Parse atom coordinates from an SDF file."""
    coords = []
    with open(lig_sdf) as f:
        lines = f.readlines()
    # SDF atom block: 3 header lines, then atom lines, then "$$$$" or "M  END"
    start = False
    for line in lines:
        if "V2000" in line or "V3000" in line:
            start = True
            continue
        if not start:
            continue
        if "$$$$" in line or "M  END" in line:
            break
        # Atom line: cols 0-9: x, 10-19: y, 20-29: z
        if len(line) >= 30 and line[0:10].strip():
            try:
                x = float(line[0:10].strip())
                y = float(line[10:20].strip())
                z = float(line[20:30].strip())
                coords.append((x, y, z))
            except ValueError:
                continue
    return coords


def compute_box(coords: list, margin: float = 10.0) -> dict:
    """Return {center_x, center_y, center_z, size_x, size_y, size_z}."""
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    return {
        "center_x": (min(xs) + max(xs)) / 2.0,
        "center_y": (min(ys) + max(ys)) / 2.0,
        "center_z": (min(zs) + max(zs)) / 2.0,
        "size_x": max((max(xs) - min(xs)) + margin * 2, 22.0),
        "size_y": max((max(ys) - min(ys)) + margin * 2, 22.0),
        "size_z": max((max(zs) - min(zs)) + margin * 2, 22.0),
    }


def get_box_for_system(sys_id: str, chain: str,
                       gt_dir: Path = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth")) -> dict | None:
    """Compute docking box for a system based on its ground-truth ligand."""
    lig_sdf = gt_dir / sys_id / "ligand_files" / f"{chain}.sdf"
    if not lig_sdf.exists():
        return None
    coords = parse_ligand_coords(lig_sdf)
    if not coords:
        return None
    return compute_box(coords, margin=10.0)
