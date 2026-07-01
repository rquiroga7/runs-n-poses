#!/bin/bash
# Check docking benchmark progress: output count per method + system count per output dir
DATASETS=/home/rquiroga/Datasets/runs-n-poses-datasets
METHOD_DIRS=(
    "$DATASETS/autodock_vina_8"
    "$DATASETS/autodock_vina_32"
    "$DATASETS/qvina_w"
    "$DATASETS/quickvina2"
    "$DATASETS/autodock_gpu"
    "$DATASETS/vina_8_meeko"
    "$DATASETS/vina_32_meeko"
    "$DATASETS/vinardo_outputs"
    "$DATASETS/vinardock_meeko"
    "$DATASETS/rdock"
    "$DATASETS/gnina"
)
METHOD_NAMES=(
    "autodock_vina_8" "autodock_vina_32" "qvina_w" "quickvina2"
    "autodock_gpu"
    "vina_8_meeko" "vina_32_meeko" "vinardock_2vinardo" "vinardock_meeko"
    "rdock" "gnina"
)

echo "=== Progress: $(date) ==="
for i in "${!METHOD_DIRS[@]}"; do
    d="${METHOD_DIRS[$i]}"
    m="${METHOD_NAMES[$i]}"
    dock_str="---"
    ana_str="---"
    if [ -d "$d" ]; then
        n=$(find "$d" \( -name "out.pdbqt" -o -name "log.csv" -o -name "output.dlg" -o -name "out.sd" -o -name "out.sdf" \) 2>/dev/null | wc -l)
        dock_str="$n"
    fi
    csv="$DATASETS/predictions/${m}.csv"
    if [ -f "$csv" ]; then
        # Count unique system_ids matching this method (handles aliases)
        n=$(python3 -c "
import csv
method='$m'
aliases={'autodock_vina_8':'vina'}
s=set()
aliases_to_check={method}
if method in aliases:
    aliases_to_check.add(aliases[method])
with open('$csv') as f:
    reader=csv.DictReader(f)
    for row in reader:
        if row.get('method','') in aliases_to_check:
            s.add(row.get('system_id','') or row.get('target',''))
print(len(s))
" 2>/dev/null)
        ana_str="$n"
    fi
    printf "%-20s %6s docked  %6s analyzed\n" "$m" "$dock_str" "$ana_str"
done
echo "=== Meeko prep ==="
for d in meeko_receptors_pdbqt meeko_receptors_pdbt meeko_ligands_pdbqt meeko_ligands_pdbt; do
    n=$(ls "$DATASETS/$d" 2>/dev/null | wc -l)
    printf "%-25s %4d systems\n" "$d" "$n"
done

echo ""
echo "=== System processes ==="
ps aux | grep -E "run_full_benchmark|03_run_vina|run_qvina|run_quickvina|run_vina_gpu|run_vina_cuda|run_autodock_gpu|run_vinardo|vinardock|mk_prepare" | grep -v grep | awk '{print $11, $12, $13}' | sort -u | head -5
