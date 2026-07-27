# PaScaL_TDMA 2.1 Study workflow

This directory contains corresponding CUDA Fortran and CUDA C++ profiling
drivers, a multi-case sweep script, curated result data, generated tables and
figures, and the current engineering report.

The two drivers use the same command-line problem definition and CSV columns.
That structural correspondence does not by itself establish numerical or
performance equivalence between the implementations. The current curated data
contains CUDA C++ rows only.

## Files

```text
../../examples/example_fortran_profile.f90      CUDA Fortran profile driver
../../examples/example_cuda_cxx_profile.cu      CUDA C++ profile driver
../../tool/run_study_sweep.sh                    case generation and execution
report_cuda_cxx_porting_tdma_solver_study.md
result/                          curated CSV data, environment, tables, figures
profile_results/tdma_profile_260702_100123.csv
                                 component-level raw CUDA C++ profile
```

The component-level raw profile is retained for traceability but is not part
of the 250-row merged dataset summarized by the curated report.

## Current evidence boundary

The checked-in merged dataset contains:

- 250 timing rows from 25 CUDA C++ cases;
- 25 CUDA C++ correctness rows;
- device-direct and host-staging MPI cases;
- rank counts 1, 2, 4, and 8;
- no CUDA Fortran timing or correctness rows.

Correctness is written only after iteration 0. The largest recorded
iteration-0 absolute error against the expected value of one is `8.860e-12`.

Both profiling drivers initialize `A`, `B`, `C`, and `D` once and then repeat
only the solver call. The solver is in-place: a one-rank solve changes `C` and
`D`, while a multi-rank solve changes `A`, `C`, and `D`. Therefore:

- iteration 0 is the only solve performed on the original initialized inputs;
- iterations 1–9 operate on successively modified inputs;
- the difference between iteration 0 and later iterations cannot be attributed
  solely to GPU warm-up or first-call overhead;
- iterations 1–9 do not measure repeated independent solves of the same
  original system.

The checked-in performance tables summarize iterations 1–9 under this current
repetition protocol. Use them only with that boundary. To measure repeated
solves of an identical mathematical problem, a future driver revision must
restore the original inputs before every measured solve.

The checked-in table headers and SVG titles use the explicit phrases
`iterations 1-9 mean` and `first call`. The current
`result/analyze_study_results.py` generator still contains the older output
labels `stable` and `warm-up`; rerunning it will overwrite those corrected
labels even though the numeric values are unchanged. The generator is left
unchanged in this documentation-only revision. Relabel regenerated assets with
the terminology above until the analysis tool is updated in a separate code
change.

## Requirements

The Study build combines both component toolchains:

- NVIDIA HPC SDK CUDA Fortran compiler and MPI Fortran wrapper;
- CUDA Toolkit with `nvcc` and an MPI C++ wrapper;
- CUDA-aware MPI for the Fortran driver and CUDA C++ device mode;
- CUDA-capable GPUs and GNU Make.

CUDA C++ host mode can run with a non-CUDA-aware MPI implementation. The
Fortran driver has no host-staging fallback.

## Build

From the repository root:

```bash
make study CUDA_ARCH=90 FC=mpifort MPICXX=mpicxx
```

This builds both implementation libraries and creates:

```text
run/example_fortran_profile
run/example_cuda_cxx_profile
```

`CUDA_ARCH=90` is the setting used for the H200 validation system, not a
universal requirement.

## Run one case

From the repository root:

```bash
mpirun -np 4 ./run/example_cuda_cxx_profile 64 64 2048 10 128 128
mpirun -np 4 ./run/example_fortran_profile 64 64 2048 10 128 128
```

Arguments for both drivers are:

```text
[n1] [n2] [n3] [iterations] [tdma_threads] [reduced_threads]
```

For the CUDA C++ driver on non-CUDA-aware MPI:

```bash
PASCAL_TDMA_MPI_MODE=host \
  mpirun -np 4 ./run/example_cuda_cxx_profile 64 64 2048 10
```

## Sweep presets

Run the script from the repository root. With no preset, it uses the
`portfolio` matrix.

Small validation subset:

```bash
STUDY_PRESET=quick ./tool/run_study_sweep.sh
```

Explicit rank/size matrix:

```bash
STUDY_PRESET=custom \
NP_LIST="1 2 4" \
SIZE_LIST="64,64,2048 128,128,4096" \
ITERATIONS=10 \
./tool/run_study_sweep.sh
```

Preview the generated manifest and launch commands without executing MPI jobs:

```bash
DRY_RUN=1 STUDY_PRESET=quick ./tool/run_study_sweep.sh
```

CUDA C++ only, host staging only:

```bash
RUN_FORTRAN=0 \
MPI_MODE=host \
MPI_MODE_LIST=host \
STUDY_PRESET=quick \
./tool/run_study_sweep.sh
```

Setting only `MPI_MODE=host` is insufficient for a host-only quick or portfolio
run because the dedicated MPI comparison cases use `MPI_MODE_LIST`, whose
default contains both `device` and `host`.

## Important sweep variables

| Variable | Default | Meaning |
| --- | --- | --- |
| `STUDY_PRESET` | `portfolio` | `portfolio`, `quick`, or `custom` case generator |
| `ITERATIONS` | `10` | solver calls per execution case |
| `RUN_FORTRAN` | `1` | include the Fortran driver |
| `RUN_CXX` | `1` | include the CUDA C++ driver |
| `CXX_DEFAULT_MPI_MODES` | value of `MPI_MODE` | C++ modes for ordinary cases |
| `MPI_MODE` | `device` | default CUDA C++ mode |
| `MPI_MODE_LIST` | `device host` | modes for MPI comparison cases |
| `NP_LIST` | `2 4 8` | rank list for `custom` |
| `SIZE_LIST` | `128,128,4096` | size list for `custom` |
| `CUDA_VISIBLE_DEVICES` | inherited | GPUs visible to each launch |
| `MPIRUN` | `mpirun` | MPI launcher command |
| `DRY_RUN` | `0` | write metadata and print commands without running |

`TDMA_THREADS`, `REDUCED_THREADS`, output paths, timestamp, scaling ranks, and
baseline rank count are also configurable near the top of
`tool/run_study_sweep.sh`.

## Generated outputs

By default, a new run writes timestamped files under `run/study/`:

```text
tdma_total_profile_<timestamp>.csv
tdma_correctness_<timestamp>.csv
tdma_case_manifest_<timestamp>.csv
tdma_environment_<timestamp>.txt
```

The manifest describes intended study cases, while the timing and correctness
files record launches that actually produced output. Always compare actual CSV
coverage with the manifest before making completeness claims. Some MPI
launchers can inherit or consume the case loop's standard input; if a sweep
stops early, use the launcher's no-stdin facility or run the printed commands
individually, then verify row and case counts.

The curated report inputs and generated assets are under
`doc/study/result/`. They are not automatically replaced by a new sweep.

## Interpreting metrics

- `total_s_max`: maximum elapsed solver time across ranks;
- `total_s_avg`: rank-average elapsed solver time;
- `compute_s_max`: maximum combined compute-phase time;
- `communication_s_max`: maximum combined MPI time;
- `packing_s_max`: maximum combined pack/unpack time;
- `throughput_Mcells_s = n1 * n2 * n3 / total_s_max / 10^6`;
- local `nrow` is either `floor(n3 / nranks)` or `ceil(n3 / nranks)`.

For cases where `n3` is divisible by the rank count, both local bounds equal
`n3 / nranks`. The report uses the maximum rank time because synchronous job
completion is limited by the slowest rank.

## Curated report

Read
[`report_cuda_cxx_porting_tdma_solver_study.md`](report_cuda_cxx_porting_tdma_solver_study.md)
for the checked-in environment and numerical tables. The report explicitly
separates iteration-0 correctness from the iterations 1–9 execution-path
timing. It does not claim a Fortran-versus-CUDA-C++ performance comparison.

## Cleanup

From the repository root:

```bash
make clean
make veryclean
```

`make clean` removes the two Study executables while preserving generated CSV
files and the curated `result/` directory. `make veryclean` additionally
removes empty generated `include/` and `lib/` directories. The repository helper
`tool/clean_for_sync.sh` provides a wider dry-run cleanup preview.
