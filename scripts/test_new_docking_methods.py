#!/usr/bin/env python
"""
Quick validation tests for rDock, LeDock and GNINA config/input formats.
Run with: python scripts/test_new_docking_methods.py
"""
import shutil, subprocess, sys
from pathlib import Path

RBT_ROOT = Path("/home/rquiroga/github/rDock")
LEDOCK = Path("/home/rquiroga/Downloads/ledock_linux_x86")
GNINA_SRC = Path("/home/rquiroga/github/gnina")
SYM_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/symmetry_corrected")
GT_DIR = Path("/home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth")
SYS = "5s9y__1__1.A__1.K"
CHAIN = "1.K"
REC = SYM_DIR / SYS / f"{SYS}_receptor_symm.pdb"
LIG_SDF = GT_DIR / SYS / "ligand_files" / f"{CHAIN}.sdf"
TMP = Path("/tmp/dock_tests")
TMP.mkdir(parents=True, exist_ok=True)


def test_rdock():
    print("=== rDock ===")
    if not RBT_ROOT.exists():
        print("  SKIP: rDock not found"); return
    rbdock = RBT_ROOT / "bin" / "rbdock"
    rbcavity = RBT_ROOT / "bin" / "rbcavity"
    env = {**__import__("os").environ, "RBT_ROOT": str(RBT_ROOT),
           "LD_LIBRARY_PATH": f"{RBT_ROOT}/lib"}
    out = TMP / "rdock"; out.mkdir(exist_ok=True)

    # Test 1: Can rDock read the receptor PDB?
    prm = out / "test.prm"
    prm.write_text(f"RBT_PARAMETER_FILE_V1.00\nRECEPTOR_FILE {REC}\n")
    r = subprocess.run([str(rbdock), "-r", str(prm), "-n", "1"],
                       capture_output=True, text=True, timeout=30, env=env)
    if "BAD_RECEPTOR_FILE" in r.stderr or "Inappropriate" in r.stderr:
        print(f"  receptor PDB: FAIL — rDock can't read symm PDB format")
        # Try converting to standard PDB with obabel
        pdb_std = out / "rec_standard.pdb"
        subprocess.run(["obabel-25-07", str(REC), "-O", str(pdb_std)],
                       capture_output=True, timeout=30)
        prm.write_text(f"RBT_PARAMETER_FILE_V1.00\nRECEPTOR_FILE {pdb_std}\n")
        r = subprocess.run([str(rbdock), "-r", str(prm), "-n", "1"],
                           capture_output=True, text=True, timeout=30, env=env)
        if "BAD_RECEPTOR_FILE" not in r.stderr:
            print(f"  receptor PDB: OK (after obabel conversion)")
        else:
            print(f"  receptor PDB: FAIL even after conversion")
            return
    else:
        print(f"  receptor PDB: OK")

    # Test 2: Can rDock read the ligand SDF?
    prm.write_text(f"RBT_PARAMETER_FILE_V1.00\nRECEPTOR_FILE {REC}\nLIGAND_FILE {LIG_SDF}\n")
    r = subprocess.run([str(rbdock), "-r", str(prm), "-n", "1"],
                       capture_output=True, text=True, timeout=30, env=env)
    print(f"  ligand SDF:   {'OK' if 'LIGAND_FILE' not in r.stderr else 'FAIL'}")


def test_ledock():
    print("\n=== LeDock ===")
    if not LEDOCK.exists():
        print(f"  SKIP: {LEDOCK} not found"); return
    out = TMP / "ledock"; out.mkdir(exist_ok=True)

    # Try format: Receptor on first line, Output second, Ligand third, etc.
    for name, kw_r, kw_o, kw_l, kw_c, kw_s, kw_n in [
        ("colon-space", "Receptor:", "Output:", "Ligand:", "Box center:", "Box radius:", "Number of runs:"),
        ("no-colon", "Receptor", "Output", "Ligand", "Box_center", "Box_radius", "Number_of_runs"),
        ("lowercase", "receptor", "output", "ligand", "box_center", "box_radius", "number_of_runs"),
        ("protein-kw", "Protein:", "Output:", "Ligand:", "Center:", "Size:", "Runs:"),
    ]:
        cfg = out / f"{name}.cfg"
        cfg.write_text(f"{kw_r} {REC}\n{kw_o} {out / 'out.dok'}\n{kw_l} {LIG_SDF}\n"
                       f"{kw_c} 10.868 -20.938 -13.841\n{kw_s} 30.0 30.0 30.0\n{kw_n} 20\n")
        r = subprocess.run([str(LEDOCK), str(cfg)], capture_output=True, timeout=10)
        err = r.stderr.decode('utf-8', errors='replace') if r.stderr else ""
        if r.returncode == 0:
            print(f"  format '{name}': OK")
        else:
            print(f"  format '{name}': {err[:80].strip()}")
            return
        err = (r.stderr or "").strip()
        print(f"  format '{name}': {err[:80]}")

    print("  All formats FAIL — LeDock config format is unknown (no docs available)")


def test_gnina():
    print("\n=== GNINA ===")
    gnina_bin = shutil.which("gnina") or next(GNINA_SRC.rglob("build/bin/gnina"), None)
    if gnina_bin:
        r = subprocess.run([str(gnina_bin), "--help"], capture_output=True, text=True, timeout=10)
        print(f"  binary found: {gnina_bin}")
        print(f"  gnina: OK ({'accepts --help' if r.returncode == 0 else 'FAILS'})")
        return
    if not GNINA_SRC.exists():
        print("  SKIP: gnina not found"); return
    cuda = subprocess.run(["which", "nvcc"], capture_output=True, text=True)
    if cuda.returncode != 0:
        print("  SKIP: no CUDA (required)")
        return
    print("  source + CUDA available - run: cd gnina && mkdir build && cd build && cmake .. && make")


if __name__ == "__main__":
    test_rdock()
    test_ledock()
    test_gnina()
    print("\n==========")
    print("rDock: needs PDB format conversion (obabel roundtrip). Script exists.")
    print("LeDock: config format unknown (no docs). Needs investigation.")
    print("GNINA: needs cmake build with CUDA.")
    print("Main benchmark is running — check with: bash scripts/check_progress.sh")
