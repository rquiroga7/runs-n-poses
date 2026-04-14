# Clean step 4 outputs
echo "=== Cleaning step 4 outputs ==="

# 1. Remove analysis JSON files
rm -f /home/rquiroga/github/runs-n-poses/examples/analysis/vinardock_2vinardo/*.json
echo "Removed analysis JSONs"

# 2. Remove predictions CSV
rm -f /home/rquiroga/Datasets/runs-n-poses-datasets/predictions/vinardock_2vinardo.csv
echo "Removed predictions CSV"

# 3. Remove intermediate PDB files (from PDBT conversion for ost)
find /home/rquiroga/Datasets/runs-n-poses-datasets/vinardo_outputs/ -name "*_model.pdb" -delete
echo "Removed model PDB files"

# 4. Remove intermediate PDB files from receptor CIF conversion for PoseBusters
find /home/rquiroga/Datasets/runs-n-poses-datasets/ground_truth/ -name "*.pdb" -delete
echo "Removed receptor PDB files"

# 5. Verify clean state
echo ""
echo "=== Verification ==="
ls /home/rquiroga/github/runs-n-poses/examples/analysis/vinardock_2vinardo/ 2>/dev/null | head -3 || echo "Analysis dir: empty"
ls /home/rquiroga/Datasets/runs-n-poses-datasets/predictions/vinardock_2vinardo.csv 2>/dev/null || echo "Predictions CSV: removed"

echo ""
echo "Step 4 output cleaned. Ready to re-run:"
echo "  python scripts/04_analyze_vinardo.py" (Clean step 4 outputs)