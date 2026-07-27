# PaScaL_TDMA 2.1 CUDA C++ port

This repository contains a CUDA C++17 + MPI port of the
[CUDA Fortran solver](fortran.md). It preserves the
many-system TDMA, rank-local reduction, `MPI_Alltoallv`, transformed reduced
solve, and full-row update flow while exposing a C++ interface.

See the [root README](../README.md) for publication, provenance, complete build
targets, citation, and repository-wide Study information.

## Implementation contract

- double-precision CUDA device arrays;
- one CUDA thread per independent tridiagonal system;
- system-contiguous flat storage compatible with the Fortran layout;
- device-buffer CUDA-aware MPI by default;
- optional host-staging MPI fallback;
- move-only RAII plan that owns CUDA work buffers and a duplicated MPI
  communicator;
- regular and phase-profiled solve entry points.

The flat array offset is:

```cpp
offset = sys + row * nsys;
```

For a three-dimensional field solved along `z`, the examples map
`sys = i + j * n1` and `row = k`.

## Requirements

- CUDA Toolkit with `nvcc`;
- an MPI C++ wrapper compatible with `nvcc -ccbin`;
- CUDA-capable GPU and compatible driver;
- C++17 support;
- GNU Make.

CUDA-aware MPI is required for the default path but not for the explicit host
fallback.

## Build

From the repository root:

```bash
make cuda-cxx CUDA_ARCH=90 MPICXX=mpicxx
```

Available targets:

| Target | Output |
| --- | --- |
| `make cuda-cxx-lib` | `lib/libpascal_tdma_cuda.a` and staged public header |
| `make cuda-cxx` | library, z-direction example, and profiling example |
| `make cuda-cxx-profile` | library and `run/ex_tdma_profile` |
| `make clean` | remove generated objects, libraries, headers/modules, and executables |

`CUDA_ARCH=90` targets the H200 system used for the current validation data.
Use the compute capability of the target GPU on another system. Override
`NVCC` or `MPICXX` when the compiler wrappers have different names or paths.

## Run the example

From the repository root:

```bash
mpirun -np 4 ./run/ex_tdma_zdirection
```

For a non-CUDA-aware MPI implementation:

```bash
PASCAL_TDMA_MPI_MODE=host \
  mpirun -np 4 ./run/ex_tdma_zdirection
```

Only the exact environment value `host` selects host staging. If the variable
is unset or has another value, device pointers are passed directly to MPI.

## Public API

```cpp
#include "pascal_tdma_cuda.hpp"

pascal_tdma::PascalTdmaPlan plan;
plan.create(nsys, MPI_COMM_WORLD, 128, 128);

pascal_tdma::solve(plan, d_a, d_b, d_c, d_d, nsys, nrow);

// Instrumented alternative. Phase boundaries add synchronization overhead.
pascal_tdma::SolveTimings timings;
pascal_tdma::solve_profiled(
    plan, d_a, d_b, d_c, d_d, nsys, nrow, &timings);

plan.destroy();
```

The last two `create` arguments are the local-TDMA and reduced-TDMA block sizes.
The plan records the selected MPI buffer mode when it is created. It is
non-copyable but movable.

`create` duplicates the supplied communicator. Destroy the plan explicitly, or
let its destructor run, before `MPI_Finalize`; this allows the duplicated
communicator to be released normally.

## Input and output arrays

`d_a`, `d_b`, `d_c`, and `d_d` must point to device allocations containing at
least `nsys * nrow` doubles in the documented layout. They represent the lower
diagonal, main diagonal, upper diagonal, and right-hand side. On return, `d_d`
contains the solution.

The solve uses caller arrays as in-place workspace:

| MPI ranks | Arrays modified |
| --- | --- |
| 1 | `d_c`, `d_d` |
| more than 1 | `d_a`, `d_c`, `d_d` |

`d_b` is not modified. Restore the original coefficients and right-hand side
before solving the same original system again. A plan may be reused for
compatible solves, but transformed arrays must not be mistaken for fresh
inputs.

## Size and decomposition constraints

- `nsys` passed to `solve` must equal the value used to create the plan.
- `nrow` must be positive for one rank.
- For multiple ranks, every rank must have `nrow >= 2`.
- The current internal distribution of reduced systems requires
  `nsys >= number of ranks`.
- `partition_1d` assigns balanced contiguous ranges. Local row lengths are
  `floor(global_nrow / nranks)` or `ceil(global_nrow / nranks)` and may differ
  by one.

The API throws C++ exceptions for detected argument, CUDA, and plan-state
errors. MPI errors inside an application should be handled consistently with
the application's MPI error policy.

## Profiling example

```bash
make cuda-cxx-profile CUDA_ARCH=90
mpirun -np 4 ./run/ex_tdma_profile 64 64 2048 10
```

Arguments:

```text
ex_tdma_profile [n1] [n2] [n3] [iterations] [tdma_threads] [reduced_threads]
```

The profiling example creates one plan and initializes its input arrays once,
then repeats only the in-place solver call. Iteration 0 is the solve of the
original initialized system. Later iterations operate on arrays modified by
earlier calls and are not independent solves of that original system.

For a component-only multi-case sweep, run the helper from any working
directory; it resolves the repository root from its own path:

```bash
NP_LIST="1 2 4 8" \
SIZE_LIST="64,64,2048 128,128,4096" \
ITERATIONS=10 \
./tool/run_tdma_profile_sweep.sh
```

The script writes CSV files under `run/profile_results/`. Set `MPI_MODE=host`
to select the host-staging path. For corresponding CUDA Fortran and CUDA C++
drivers, use the [Study workflow](study/README.md).

## Files

- `src/pascal_tdma_cuda.hpp`: tracked public API and plan definition;
- `include/pascal_tdma_cuda.hpp`: generated staged header;
- `src/pascal_tdma_cuda.cu`: kernels, plan implementation, and MPI wrapper;
- `examples/ex_tdma_zdirection.cu`: standalone z-direction example;
- `examples/ex_tdma_profile.cu`: component profiling example;
- `tool/run_tdma_profile_sweep.sh`: component sweep helper;
- [`cuda_cxx_porting_notes.md`](cuda_cxx_porting_notes.md): preserved porting contract and design
  decisions.

## Validation boundary

The port and current Study workflow were built and run on an 8-GPU NVIDIA H200
system with CUDA 12.9 and Open MPI. The curated dataset currently contains
CUDA C++ rows only. See the [Study report](study/report_cuda_cxx_porting_tdma_solver_study.md)
for the exact environment, case coverage, numerical results, and timing-method
caveat.

## License and citation

This component is distributed under the repository's [MIT License](../LICENSE).
Please cite the PaScaL_TDMA 2.1 paper described in the
[root citation section](../README.md#citation).
