#!/usr/bin/env bash
set -euo pipefail

MPIRUN=${MPIRUN:-mpirun}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
FORTRAN_EXE=${FORTRAN_EXE:-"$ROOT_DIR/run/example_fortran_profile"}
CXX_EXE=${CXX_EXE:-"$ROOT_DIR/run/example_cuda_cxx_profile"}

STUDY_PRESET=${STUDY_PRESET:-portfolio}
ITERATIONS=${ITERATIONS:-10}
TDMA_THREADS=${TDMA_THREADS:-128}
REDUCED_THREADS=${REDUCED_THREADS:-128}
BASELINE_NP=${BASELINE_NP:-2}
SCALING_NP_LIST=${SCALING_NP_LIST:-"2 4 8"}
RUN_NP1_REFERENCE=${RUN_NP1_REFERENCE:-1}

# C++ uses CUDA-aware MPI device buffers by default. Host mode is added only
# for the dedicated MPI-mode comparison cases unless the user overrides this.
MPI_MODE=${MPI_MODE:-device}
CXX_DEFAULT_MPI_MODES=${CXX_DEFAULT_MPI_MODES:-"$MPI_MODE"}
MPI_MODE_LIST=${MPI_MODE_LIST:-"device host"}

RUN_FORTRAN=${RUN_FORTRAN:-1}
RUN_CXX=${RUN_CXX:-1}
DRY_RUN=${DRY_RUN:-0}
TIMESTAMP=${TIMESTAMP:-$(date +%y%m%d_%H%M%S)}
OUT=${OUT:-"$ROOT_DIR/run/study/tdma_total_profile_${TIMESTAMP}.csv"}
CORRECTNESS_OUT=${CORRECTNESS_OUT:-"$ROOT_DIR/run/study/tdma_correctness_${TIMESTAMP}.csv"}
ENV_OUT=${ENV_OUT:-"$ROOT_DIR/run/study/tdma_environment_${TIMESTAMP}.txt"}
MANIFEST_OUT=${MANIFEST_OUT:-"$ROOT_DIR/run/study/tdma_case_manifest_${TIMESTAMP}.csv"}

# Custom preset runs the explicit NP_LIST x SIZE_LIST matrix supplied by the
# user.
NP_LIST=${NP_LIST:-"2 4 8"}
SIZE_LIST=${SIZE_LIST:-"128,128,4096"}

usage() {
    cat <<'USAGE'
Usage:
  ./tool/run_study_sweep.sh
  STUDY_PRESET=quick ./tool/run_study_sweep.sh
  STUDY_PRESET=custom NP_LIST="2 4 8" SIZE_LIST="128,128,4096" ./tool/run_study_sweep.sh

Study presets:
  portfolio  Full benchmark study matrix:
             correctness, phase breakdown, strong scaling, weak scaling,
             nsys/nrow sensitivity, MPI device-vs-host comparison.
  quick      Small validation subset of the benchmark matrix.
  custom     Custom NP_LIST x SIZE_LIST execution.

Important variables:
  BASELINE_NP=2
  SCALING_NP_LIST="2 4 8"
  ITERATIONS=10
  CXX_DEFAULT_MPI_MODES="device"
  MPI_MODE_LIST="device host"    # used only by mpi_mode_compare cases
  RUN_FORTRAN=1
  RUN_CXX=1
  DRY_RUN=1                      # write manifest/environment, print commands only
USAGE
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    usage
    exit 0
elif [[ $# -gt 0 ]]; then
    echo "error: unexpected argument: $1" >&2
    usage >&2
    exit 2
fi

case "$STUDY_PRESET" in
    portfolio|quick|custom)
        ;;
    *)
        echo "error: unknown STUDY_PRESET=$STUDY_PRESET" >&2
        usage >&2
        exit 2
        ;;
esac

if [[ "$RUN_FORTRAN" != "1" && "$RUN_CXX" != "1" ]]; then
    echo "error: at least one of RUN_FORTRAN or RUN_CXX must be 1" >&2
    exit 2
fi

if [[ "$DRY_RUN" != "1" ]]; then
    if [[ "$RUN_FORTRAN" == "1" && ! -x "$FORTRAN_EXE" ]]; then
        echo "error: Fortran profile driver not found: $FORTRAN_EXE (run 'make study')" >&2
        exit 1
    fi
    if [[ "$RUN_CXX" == "1" && ! -x "$CXX_EXE" ]]; then
        echo "error: CUDA C++ profile driver not found: $CXX_EXE (run 'make study')" >&2
        exit 1
    fi
fi

mkdir -p "$(dirname "$OUT")"
mkdir -p "$(dirname "$CORRECTNESS_OUT")"
mkdir -p "$(dirname "$ENV_OUT")"
mkdir -p "$(dirname "$MANIFEST_OUT")"

EXEC_CASES="$(mktemp "${TMPDIR:-/tmp}/tdma_exec_cases.XXXXXX")"
trap 'rm -f "$EXEC_CASES"' EXIT

append_timing_csv() {
    if [[ ! -s "$OUT" ]]; then
        tee -a "$OUT"
    else
        awk 'NR == 1 && /^solver,/ { next } { print }' | tee -a "$OUT"
    fi
}

manifest_header() {
    cat > "$MANIFEST_OUT" <<'EOF'
study_suite,case_id,nranks,n1,n2,n3,baseline_nranks,scaling_kind,cxx_mpi_modes,notes
EOF
}

add_case() {
    local suite="$1"
    local np="$2"
    local n1="$3"
    local n2="$4"
    local n3="$5"
    local scaling_kind="$6"
    local cxx_modes="$7"
    local notes="$8"
    local case_id="${suite}_np${np}_${n1}x${n2}x${n3}"
    local mode

    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$suite" "$case_id" "$np" "$n1" "$n2" "$n3" \
        "$BASELINE_NP" "$scaling_kind" "$cxx_modes" "$notes" >> "$MANIFEST_OUT"

    if [[ "$RUN_FORTRAN" == "1" ]]; then
        printf 'fortran-original,device,%s,%s,%s,%s\n' "$np" "$n1" "$n2" "$n3" >> "$EXEC_CASES"
    fi

    if [[ "$RUN_CXX" == "1" ]]; then
        for mode in $cxx_modes; do
            printf 'cuda-cxx,%s,%s,%s,%s,%s\n' "$mode" "$np" "$n1" "$n2" "$n3" >> "$EXEC_CASES"
        done
    fi
}

add_size_case() {
    local suite="$1"
    local np="$2"
    local size="$3"
    local scaling_kind="$4"
    local cxx_modes="$5"
    local notes="$6"
    local n1 n2 n3

    IFS=',' read -r n1 n2 n3 <<< "$size"
    add_case "$suite" "$np" "$n1" "$n2" "$n3" "$scaling_kind" "$cxx_modes" "$notes"
}

build_quick_cases() {
    local np

    if [[ "$RUN_NP1_REFERENCE" == "1" ]]; then
        add_size_case "single_gpu_reference" 1 "128,128,4096" "reference" \
            "$CXX_DEFAULT_MPI_MODES" "local_tdma_no_mpi"
    fi

    for np in 2 4; do
        add_size_case "strong_scaling" "$np" "128,128,4096" "strong_2gpu_baseline" \
            "$CXX_DEFAULT_MPI_MODES" "fixed_global_problem"
    done

    add_size_case "weak_nrow_scaling" 2 "128,128,2048" "weak_nrow_2gpu_baseline" \
        "$CXX_DEFAULT_MPI_MODES" "fixed_nsys_n3_scales_with_rank"
    add_size_case "weak_nrow_scaling" 4 "128,128,4096" "weak_nrow_2gpu_baseline" \
        "$CXX_DEFAULT_MPI_MODES" "fixed_nsys_n3_scales_with_rank"

    for np in 2 4; do
        add_size_case "mpi_mode_compare" "$np" "128,128,4096" "mpi_mode" \
            "$MPI_MODE_LIST" "cxx_device_vs_host_for_same_case"
    done
}

build_custom_cases() {
    local np size
    for np in $NP_LIST; do
        for size in $SIZE_LIST; do
            add_size_case "custom" "$np" "$size" "custom" \
                "$CXX_DEFAULT_MPI_MODES" "custom_case"
        done
    done
}

build_portfolio_cases() {
    local np

    # Correctness is collected for every execution case. These reference cases
    # keep a single-GPU local Thomas baseline in the dataset without making it
    # the scaling baseline.
    if [[ "$RUN_NP1_REFERENCE" == "1" ]]; then
        add_size_case "single_gpu_reference" 1 "64,64,2048" "reference" \
            "$CXX_DEFAULT_MPI_MODES" "local_tdma_no_mpi_small"
        add_size_case "single_gpu_reference" 1 "128,128,4096" "reference" \
            "$CXX_DEFAULT_MPI_MODES" "local_tdma_no_mpi_medium"
    fi

    # Strong scaling: fixed global problem size, baseline is np=2.
    for np in $SCALING_NP_LIST; do
        add_size_case "strong_scaling" "$np" "128,128,4096" "strong_2gpu_baseline" \
            "$CXX_DEFAULT_MPI_MODES" "fixed_global_medium"
        add_size_case "strong_scaling" "$np" "256,256,4096" "strong_2gpu_baseline" \
            "$CXX_DEFAULT_MPI_MODES" "fixed_global_nsys_rich"
    done

    # Weak scaling A: keep nsys fixed and scale n3 with rank count so nrow per
    # rank stays constant. This isolates global line-length growth.
    add_size_case "weak_nrow_scaling" 2 "128,128,2048" "weak_nrow_2gpu_baseline" \
        "$CXX_DEFAULT_MPI_MODES" "fixed_nsys_local_nrow_1024"
    add_size_case "weak_nrow_scaling" 4 "128,128,4096" "weak_nrow_2gpu_baseline" \
        "$CXX_DEFAULT_MPI_MODES" "fixed_nsys_local_nrow_1024"
    add_size_case "weak_nrow_scaling" 8 "128,128,8192" "weak_nrow_2gpu_baseline" \
        "$CXX_DEFAULT_MPI_MODES" "fixed_nsys_local_nrow_1024"

    # Weak scaling B: keep n3 fixed and scale nsys with rank count so local
    # nsys*nrow work stays roughly constant. This isolates system-count growth.
    add_size_case "weak_nsys_scaling" 2 "128,128,2048" "weak_nsys_2gpu_baseline" \
        "$CXX_DEFAULT_MPI_MODES" "n2_scaled_local_work_constant"
    add_size_case "weak_nsys_scaling" 4 "128,256,2048" "weak_nsys_2gpu_baseline" \
        "$CXX_DEFAULT_MPI_MODES" "n2_scaled_local_work_constant"
    add_size_case "weak_nsys_scaling" 8 "128,512,2048" "weak_nsys_2gpu_baseline" \
        "$CXX_DEFAULT_MPI_MODES" "n2_scaled_local_work_constant"

    # Sensitivity studies around the 2-GPU baseline and the largest rank count.
    for np in 2 8; do
        add_size_case "nsys_sensitivity" "$np" "64,64,4096" "nsys_sweep" \
            "$CXX_DEFAULT_MPI_MODES" "n3_fixed_vary_nsys"
        add_size_case "nsys_sensitivity" "$np" "128,128,4096" "nsys_sweep" \
            "$CXX_DEFAULT_MPI_MODES" "n3_fixed_vary_nsys"
        add_size_case "nsys_sensitivity" "$np" "128,256,4096" "nsys_sweep" \
            "$CXX_DEFAULT_MPI_MODES" "n3_fixed_vary_nsys"

        add_size_case "nrow_sensitivity" "$np" "128,128,2048" "nrow_sweep" \
            "$CXX_DEFAULT_MPI_MODES" "nsys_fixed_vary_n3"
        add_size_case "nrow_sensitivity" "$np" "128,128,4096" "nrow_sweep" \
            "$CXX_DEFAULT_MPI_MODES" "nsys_fixed_vary_n3"
        add_size_case "nrow_sensitivity" "$np" "128,128,8192" "nrow_sweep" \
            "$CXX_DEFAULT_MPI_MODES" "nsys_fixed_vary_n3"
    done

    # MPI mode comparison: C++ device-buffer path vs host-staging fallback.
    for np in $SCALING_NP_LIST; do
        add_size_case "mpi_mode_compare" "$np" "128,128,4096" "mpi_mode" \
            "$MPI_MODE_LIST" "cxx_device_vs_host_for_same_case"
    done
}

capture_environment() {
    {
        echo "# PaScaL_TDMA 2.1 Study Environment"
        echo "date=$(date '+%Y-%m-%dT%H:%M:%S%z')"
        echo "hostname=$(hostname)"
        echo "pwd=$PWD"
        echo "root_dir=$ROOT_DIR"
        echo "script_dir=$SCRIPT_DIR"
        echo "study_preset=$STUDY_PRESET"
        echo "baseline_np=$BASELINE_NP"
        echo "scaling_np_list=$SCALING_NP_LIST"
        echo "run_np1_reference=$RUN_NP1_REFERENCE"
        echo "custom_np_list=$NP_LIST"
        echo "custom_size_list=$SIZE_LIST"
        echo "iterations=$ITERATIONS"
        echo "tdma_threads=$TDMA_THREADS"
        echo "reduced_threads=$REDUCED_THREADS"
        echo "cxx_default_mpi_modes=$CXX_DEFAULT_MPI_MODES"
        echo "mpi_mode_list=$MPI_MODE_LIST"
        echo "run_fortran=$RUN_FORTRAN"
        echo "run_cxx=$RUN_CXX"
        echo "dry_run=$DRY_RUN"
        echo "cuda_visible_devices=${CUDA_VISIBLE_DEVICES:-}"
        echo "timing_csv=$OUT"
        echo "correctness_csv=$CORRECTNESS_OUT"
        echo "environment_file=$ENV_OUT"
        echo "case_manifest=$MANIFEST_OUT"
        echo
        echo "## git"
        git -C "$ROOT_DIR" rev-parse HEAD 2>/dev/null || true
        git -C "$ROOT_DIR" status --short 2>/dev/null || true
        echo
        echo "## nvidia-smi"
        if command -v nvidia-smi >/dev/null 2>&1; then
            nvidia-smi || true
            echo
            echo "## nvidia-smi topo -m"
            nvidia-smi topo -m || true
        else
            echo "nvidia-smi not found"
        fi
        echo
        echo "## nvcc --version"
        if command -v nvcc >/dev/null 2>&1; then
            nvcc --version || true
        else
            echo "nvcc not found"
        fi
        echo
        echo "## mpirun --version"
        "$MPIRUN" --version || true
        echo
        echo "## mpifort --version"
        if command -v mpifort >/dev/null 2>&1; then
            mpifort --version || true
        else
            echo "mpifort not found"
        fi
        echo
        echo "## mpicxx --version"
        if command -v mpicxx >/dev/null 2>&1; then
            mpicxx --version || true
        else
            echo "mpicxx not found"
        fi
    } > "$ENV_OUT"
}

manifest_header
case "$STUDY_PRESET" in
    portfolio)
        build_portfolio_cases
        ;;
    quick)
        build_quick_cases
        ;;
    custom)
        build_custom_cases
        ;;
    *)
        echo "error: unknown STUDY_PRESET=$STUDY_PRESET" >&2
        usage >&2
        exit 2
        ;;
esac

capture_environment

sort -u "$EXEC_CASES" | while IFS=, read -r implementation mode np n1 n2 n3; do
    echo "running implementation=$implementation mpi_mode=$mode np=$np n1=$n1 n2=$n2 n3=$n3 iterations=$ITERATIONS" >&2

    if [[ "$implementation" == "fortran-original" ]]; then
        if [[ "$DRY_RUN" == "1" ]]; then
            echo "[dry-run] $MPIRUN -np $np $FORTRAN_EXE $n1 $n2 $n3 $ITERATIONS $TDMA_THREADS $REDUCED_THREADS" >&2
            continue
        fi
        PASCAL_TDMA_CORRECTNESS_OUT="$CORRECTNESS_OUT" \
        "$MPIRUN" -np "$np" "$FORTRAN_EXE" "$n1" "$n2" "$n3" \
            "$ITERATIONS" "$TDMA_THREADS" "$REDUCED_THREADS" < /dev/null \
            | append_timing_csv
    elif [[ "$implementation" == "cuda-cxx" ]]; then
        if [[ "$DRY_RUN" == "1" ]]; then
            echo "[dry-run] PASCAL_TDMA_MPI_MODE=$mode $MPIRUN -np $np $CXX_EXE $n1 $n2 $n3 $ITERATIONS $TDMA_THREADS $REDUCED_THREADS" >&2
            continue
        fi
        PASCAL_TDMA_CORRECTNESS_OUT="$CORRECTNESS_OUT" \
        PASCAL_TDMA_MPI_MODE="$mode" \
        "$MPIRUN" -np "$np" "$CXX_EXE" "$n1" "$n2" "$n3" \
            "$ITERATIONS" "$TDMA_THREADS" "$REDUCED_THREADS" < /dev/null \
            | append_timing_csv
    else
        echo "error: unknown implementation=$implementation" >&2
        exit 2
    fi
done

if [[ "$DRY_RUN" == "1" ]]; then
    echo "dry-run: timing and correctness CSV files were not generated" >&2
else
    echo "wrote $OUT" >&2
    echo "wrote $CORRECTNESS_OUT" >&2
fi
echo "wrote $ENV_OUT" >&2
echo "wrote $MANIFEST_OUT" >&2
