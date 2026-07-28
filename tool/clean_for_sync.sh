#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
DRY_RUN=1

usage() {
    cat <<'USAGE'
Usage:
  tool/clean_for_sync.sh          # preview generated files
  tool/clean_for_sync.sh --apply  # remove the previewed files

This helper removes build products and transient profiler output while
preserving source files, scripts, CSV/TXT study data, and doc/study/result/.
The regular `make clean` target is narrower and is preferred after builds.
USAGE
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
    usage
    exit 0
elif [[ "${1:-}" == "--apply" ]]; then
    DRY_RUN=0
elif [[ $# -gt 0 ]]; then
    usage >&2
    exit 2
fi

remove_path() {
    local path="$1"
    [[ -e "$path" ]] || return 0

    if [[ "$DRY_RUN" -eq 1 ]]; then
        printf '[dry-run] remove %s\n' "${path#$ROOT_DIR/}"
    else
        rm -rf -- "$path"
        printf 'removed %s\n' "${path#$ROOT_DIR/}"
    fi
}

# Exact products of the checked-in Makefiles.
remove_path "$ROOT_DIR/src/obj"
remove_path "$ROOT_DIR/lib/libpascal_tdma.a"
remove_path "$ROOT_DIR/lib/libpascal_tdma_cuda.a"
remove_path "$ROOT_DIR/include/pascal_tdma_cuda.hpp"
remove_path "$ROOT_DIR/run/a.out"
remove_path "$ROOT_DIR/run/ex_tdma_zdirection"
remove_path "$ROOT_DIR/run/ex_tdma_profile"
remove_path "$ROOT_DIR/run/example_fortran_profile"
remove_path "$ROOT_DIR/run/example_cuda_cxx_profile"

if [[ -d "$ROOT_DIR/include" ]]; then
    while IFS= read -r -d '' path; do
        remove_path "$path"
    done < <(find "$ROOT_DIR/include" -maxdepth 1 -type f \
        \( -name '*.mod' -o -name '*.smod' \) -print0)
fi

# Wider transient-file cleanup. Curated documentation and all CSV/TXT data are
# excluded deliberately; the helper is safe to preview before applying.
while IFS= read -r -d '' path; do
    remove_path "$path"
done < <(
    find "$ROOT_DIR" \
        -path "$ROOT_DIR/.git" -prune -o \
        -path "$ROOT_DIR/src/obj" -prune -o \
        -path "$ROOT_DIR/include" -prune -o \
        -path "$ROOT_DIR/doc/study/result" -prune -o \
        -type f \( \
            -name '*.o' -o \
            -name '*.mod' -o \
            -name '*.smod' -o \
            -name '*.so' -o \
            -name '*.dylib' -o \
            -name '*.nsys-rep' -o \
            -name '*.qdrep' -o \
            -name '*.ncu-rep' -o \
            -name '*.sqlite' -o \
            -name 'log_*.out' -o \
            -name '*.log' -o \
            -name '*.err' -o \
            -name '*.vtk' -o \
            -name '*.vtr' -o \
            -name '*.vtu' -o \
            -name '*.pvtr' -o \
            -name '*.pvtu' -o \
            -name '*.pvd' -o \
            -name '*.dat' -o \
            -name '*.bin' -o \
            -name '*.plt' -o \
            -name '.DS_Store' \
        \) ! -name '*.csv' ! -name '*.txt' -print0
)

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "dry-run only. Re-run with --apply after reviewing the list."
fi
