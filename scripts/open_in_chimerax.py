#!/usr/bin/env python3
"""Open ground-truth and docking outputs for a system in ChimeraX.

Usage examples:
  python scripts/open_in_chimerax.py 7ftv --dry-run
  python scripts/open_in_chimerax.py 7ftv --chimerax-path /usr/bin/chimerax

The script looks up the full `group_key` in `scripts/single_ligand_systems.txt`,
collects receptor/ligand and docking output files under `runs-n-poses-datasets`,
and launches ChimeraX with those files. Use `--dry-run` to only print paths.
"""
from pathlib import Path
import argparse
import subprocess
import sys


def find_group_key(code: str, systems_file: Path):
    code = code.strip()
    if not systems_file.exists():
        return None
    with systems_file.open() as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            # lines typically start with the 4-letter code
            if s.startswith(code) or s.split('__', 1)[0] == code:
                return s
    return None


def find_files(repo_root: Path, group_key: str):
    dataset = repo_root / 'runs-n-poses-datasets'
    found = []

    # ground truth
    gt = dataset / 'ground_truth' / group_key
    if gt.exists():
        for name in ('receptor.pdb', 'receptor.cif', 'receptor.pdbqt'):
            p = gt / name
            if p.exists():
                found.append(p)
        ligdir = gt / 'ligand_files'
        if ligdir.exists():
            for p in sorted(ligdir.glob('*')):
                if p.is_file():
                    found.append(p)

    # docking output locations to check (common names in this repo)
    candidate_roots = [
        'autodock_vina_32', 'autodock_vina_8', 'autodock_vina',
        'vinardo_outputs', 'vinardock_2vinardo', 'vinardo'
    ]
    for root_name in candidate_roots:
        base = dataset / root_name / group_key
        if base.exists():
            # collect typical file types
            for ext in ('*.pdbqt', '*.pdb', '*.mol2', '*.sdf'):
                for p in sorted(base.rglob(ext)):
                    if p.is_file():
                        found.append(p)

    # fallback: search entire dataset for directories named group_key
    if not found and dataset.exists():
        for p in dataset.rglob(group_key):
            if p.is_dir():
                for child in sorted(p.rglob('*')):
                    if child.is_file() and child.suffix.lower() in ('.pdbqt', '.pdb', '.sdf', '.mol2'):
                        found.append(child)

    # preserve order and unique
    seen = set()
    uniq = []
    for p in found:
        s = str(p)
        if s not in seen:
            seen.add(s)
            uniq.append(p)
    return uniq


def main():
    p = argparse.ArgumentParser()
    p.add_argument('code', help='4-letter PDB-like code to look up (e.g. 7ftv)')
    p.add_argument('--dry-run', action='store_true', help='Only print found files and the ChimeraX command')
    p.add_argument('--chimerax-path', default='chimerax', help='ChimeraX executable or path')
    p.add_argument('--systems-file', default=str(Path(__file__).resolve().parents[1] / 'scripts' / 'single_ligand_systems.txt'))
    args = p.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    systems_file = Path(args.systems_file)

    group_key = find_group_key(args.code, systems_file)
    if not group_key:
        print(f'ERROR: code "{args.code}" not found in {systems_file}', file=sys.stderr)
        sys.exit(2)

    files = find_files(repo_root, group_key)

    if not files:
        print(f'No files found for group_key {group_key} under runs-n-poses-datasets', file=sys.stderr)
        sys.exit(3)

    print('Found group_key:', group_key)
    print('Files to open (in order):')
    for f in files:
        print('  ', f)

    cmd = [args.chimerax_path] + [str(x) for x in files]

    if args.dry_run:
        print('\nChimeraX command:')
        print(' '.join(cmd))
        return

    # Launch ChimeraX with the files
    try:
        print('Launching ChimeraX...')
        subprocess.run(cmd)
    except FileNotFoundError:
        print('ERROR: ChimeraX executable not found at', args.chimerax_path, file=sys.stderr)
        print('Use --chimerax-path to point to your local ChimeraX binary.', file=sys.stderr)
        sys.exit(4)


if __name__ == '__main__':
    main()
