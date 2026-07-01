#!/bin/bash
# Background runner for the docking benchmark. Logs go to logs/benchmark/.
set -e

cd /home/rquiroga/github/runs-n-poses

# Set number of CPU threads
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8

# Use the venv Python for Meeko calls
export PATH="/home/rquiroga/github/runs-n-poses/.venv/bin:$PATH"

LOG_DIR=/home/rquiroga/github/runs-n-poses/logs/benchmark
mkdir -p "$LOG_DIR"

# Run the full benchmark in the background
nohup /home/rquiroga/github/runs-n-poses/.venv/bin/python scripts/run_full_benchmark.py \
    --steps receptor,ligand,dock,analyze \
    --threads 8 \
    --resume \
    --log-dir "$LOG_DIR" \
    > "$LOG_DIR/nohup.out" 2>&1 &
BG_PID=$!

echo "Started benchmark in background as PID $BG_PID"
echo "Logs: $LOG_DIR/benchmark.log"
echo "stdout: $LOG_DIR/nohup.out"
echo "$BG_PID" > "$LOG_DIR/pid"

# Tail the log to show progress for a few seconds
sleep 10
tail -30 "$LOG_DIR/benchmark.log" 2>/dev/null || echo "Log not yet created"
