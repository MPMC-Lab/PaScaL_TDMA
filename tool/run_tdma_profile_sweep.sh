#!/usr/bin/env bash
set -euo pipefail

MPIRUN=${MPIRUN:-mpirun}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
EXE=${EXE:-"$ROOT_DIR/run/ex_tdma_profile"}
NP_LIST=${NP_LIST:-"1 2 4 8"}
SIZE_LIST=${SIZE_LIST:-"64,64,2048 128,128,2048 128,128,4096"}
ITERATIONS=${ITERATIONS:-10}
TDMA_THREADS=${TDMA_THREADS:-128}
REDUCED_THREADS=${REDUCED_THREADS:-128}
MPI_MODE=${MPI_MODE:-device}
OUT=${OUT:-"$ROOT_DIR/run/profile_results/tdma_profile_$(date +%y%m%d_%H%M%S).csv"}

if [[ ! -x "$EXE" ]]; then
    echo "error: executable not found: $EXE (run 'make cuda-cxx-profile')" >&2
    exit 1
fi

mkdir -p "$(dirname "$OUT")"

append_csv() {
    if [[ ! -s "$OUT" ]]; then
        tee -a "$OUT"
    else
        awk 'NR == 1 && /^solver,/ { next } { print }' | tee -a "$OUT"
    fi
}

for np in $NP_LIST; do
    for size in $SIZE_LIST; do
        IFS=',' read -r n1 n2 n3 <<< "$size"
        echo "running np=$np n1=$n1 n2=$n2 n3=$n3 iterations=$ITERATIONS mode=$MPI_MODE" >&2

        if [[ "$MPI_MODE" == "default" ]]; then
            "$MPIRUN" -np "$np" "$EXE" "$n1" "$n2" "$n3" "$ITERATIONS" "$TDMA_THREADS" "$REDUCED_THREADS" \
                | append_csv
        else
            PASCAL_TDMA_MPI_MODE="$MPI_MODE" \
            "$MPIRUN" -np "$np" "$EXE" "$n1" "$n2" "$n3" "$ITERATIONS" "$TDMA_THREADS" "$REDUCED_THREADS" \
                | append_csv
        fi
    done
done

echo "wrote $OUT" >&2
