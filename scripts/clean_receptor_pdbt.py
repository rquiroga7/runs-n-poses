#!/usr/bin/env python
"""
Post-process receptor PDBT files to remove unsupported atom types.

Vinardock crashes on receptors containing certain metal/element types.
This script strips those atoms from the PDBT files.
"""

import argparse
from pathlib import Path

# Unsupported atom types that cause vinardock crashes
UNSUPPORTED = {
    'V', 'MO', 'SE', 'B', 'I', 'AT', 'U', 'PU', 'AM', 'CM', 'BK', 'CF',
    'ES', 'FM', 'MD', 'NO', 'LR', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT',
    'DS', 'RG', 'CN', 'FL', 'LV', 'TS', 'OG', 'HG', 'Hg', 'Ni', 'Na'
}


def clean_pdbt(input_file: str, output_file: str = None) -> dict:
    """
    Remove lines with unsupported atom types from a PDBT file.
    Returns stats dict.
    """
    if output_file is None:
        output_file = input_file

    total = 0
    removed = 0
    removed_types = set()
    kept_lines = []

    with open(input_file) as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            total += 1
            atype = line[76:].strip()
            if atype in UNSUPPORTED:
                removed += 1
                removed_types.add(atype)
                continue
        kept_lines.append(line)

    with open(output_file, 'w') as f:
        f.writelines(kept_lines)

    return {
        "total": total,
        "removed": removed,
        "removed_types": removed_types,
        "kept": total - removed,
    }


def main():
    parser = argparse.ArgumentParser(description="Clean unsupported atoms from receptor PDBT files")
    parser.add_argument(
        "--receptor-dir",
        type=str,
        default="/home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_inputs/receptors",
        help="Directory containing receptor PDBT files"
    )
    parser.add_argument(
        "--inplace",
        action="store_true",
        help="Modify files in-place (default: dry run)"
    )
    parser.add_argument(
        "--system-id",
        type=str,
        default=None,
        help="Process only a specific system ID"
    )

    args = parser.parse_args()

    receptor_dir = Path(args.receptor_dir)

    if args.system_id:
        pdbt_files = list(receptor_dir.glob(f"{args.system_id}/*receptor*.pdbt"))
    else:
        pdbt_files = list(receptor_dir.glob("**/*receptor*.pdbt"))

    print(f"Found {len(pdbt_files)} receptor PDBT files")

    total_removed = 0
    systems_cleaned = 0

    for pdbt in pdbt_files:
        stats = clean_pdbt(str(pdbt))
        if stats["removed"] > 0:
            if args.inplace:
                # Already written in clean_pdbt when output_file == input_file
                systems_cleaned += 1
                total_removed += stats["removed"]
                print(f"  Cleaned {pdbt.name}: removed {stats['removed']} atoms ({stats['removed_types']})")
            else:
                print(f"  {pdbt.name}: would remove {stats['removed']} atoms ({stats['removed_types']})")

    if not args.inplace:
        print("\nDry run. Use --inplace to actually modify files.")
    else:
        print(f"\nDone! Removed {total_removed} unsupported atoms from {systems_cleaned} files")


if __name__ == "__main__":
    main()
