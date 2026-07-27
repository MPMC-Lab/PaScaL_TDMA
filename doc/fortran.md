# PaScaL_TDMA 2.1 CUDA Fortran implementation

This manual describes the CUDA Fortran + MPI implementation of the PaScaL_TDMA
2.1 multi-GPU tridiagonal solver. The maintained source includes changes made
after the version 3 publication archive, including a nonblocking-request fix
and optional phase-level profiling support.

For the paper, archive, repository-wide provenance, and citation information,
see the [root README](../README.md).

## Features

- solves many independent double-precision tridiagonal systems on GPUs;
- decomposes the row direction across MPI ranks;
- reduces each rank-local segment to two interface rows;
- redistributes reduced systems with CUDA-aware `MPI_Alltoallv`;
- exposes reusable plan, solve, profiled-solve, and cleanup procedures;
- includes a three-dimensional z-direction example.

The Fortran implementation always passes device buffers to MPI. It therefore
requires a CUDA-aware MPI implementation and does not provide a host-staging
fallback.

## Directory layout

```text
Makefile                    repository build targets
Makefile.inc                checked-in compiler and CUDA flag defaults
src/pascal_tdma_cuda.f90    module and solver implementation
examples/ex_tdma_zdirection.f90
lib/libpascal_tdma.a        generated static library
include/*.mod               generated Fortran module files
run/a.out                   generated example executable
```

## Requirements

- NVIDIA HPC SDK with CUDA Fortran support;
- an MPI Fortran wrapper that invokes the NVIDIA compiler;
- CUDA-aware MPI;
- a CUDA-capable GPU and compatible driver;
- GNU Make.

The checked-in `Makefile.inc` defaults to `FC=mpifort` and `CUDA_ARCH=90`.
`FC` must wrap NVIDIA `nvfortran`. Override the wrapper and architecture on the
command line for the target system.

## Build

From the repository root:

```bash
make all FC=mpifort CUDA_ARCH=90
```

Available targets:

| Target | Output |
| --- | --- |
| `make lib` | `lib/libpascal_tdma.a` and `include/*.mod` |
| `make example` | library plus `run/a.out` |
| `make all` | same component outputs as `make example` |
| `make clean` | remove generated objects, modules, library contents, and executables |
| `make veryclean` | also remove empty generated `include/` and `lib/` directories |

`make clean` retains the `lib/`, `include/`, and `run/` directories. Both
cleanup targets preserve `run/job.sh`.

The explicit implementation alias is:

```bash
make fortran FC=mpifort CUDA_ARCH=90
```

## Run the example

From the repository root, after building:

```bash
mpirun -np 4 ./run/a.out
```

The example creates `64 x 64` independent systems with global row length
`2048`, partitions the row direction across ranks, assigns each rank to
`rank % visible_device_count`, and solves a problem whose expected solution is
one.

`run/job.sh` is a site-specific Slurm example. Submit it from the generated
run directory so `SLURM_SUBMIT_DIR` contains `a.out`:

```bash
cd run
sbatch job.sh
```

Avoid accidental GPU oversubscription unless the MPI placement and application
are configured for it.

## Public solver interface

```fortran
use PaScaL_TDMA_cuda

type(ptdma_plan_cuda)   :: plan
type(ptdma_timing_cuda) :: timing

call pascal_plan_create(plan, nsys, comm, rank, nranks, &
                        tdma_threads, reduced_threads)

call pascal_solver(plan, A_d, B_d, C_d, D_d, nsys, nrow)

! Instrumented alternative; it follows the same numerical path but adds
! synchronization around measured phases.
call pascal_solver_profiled(plan, A_d, B_d, C_d, D_d, &
                            nsys, nrow, timing)

call pascal_plan_clean(plan)
```

Call `pascal_plan_clean` after the last solve and before `MPI_Finalize`.
`pascal_setcudathread` exists in the module but is currently an empty
placeholder; configure the two supported block sizes through
`pascal_plan_create`.

## Array layout and input contract

The solver arguments are double-precision CUDA Fortran device arrays:

```fortran
real(8), device :: A_d(0:nsys-1,0:nrow-1)
real(8), device :: B_d(0:nsys-1,0:nrow-1)
real(8), device :: C_d(0:nsys-1,0:nrow-1)
real(8), device :: D_d(0:nsys-1,0:nrow-1)
```

The system index is contiguous in memory. `A`, `B`, and `C` are the lower,
main, and upper diagonals; `D` is the right-hand side on entry and the solution
on return.

The solve is in-place:

| MPI ranks | Arrays modified |
| --- | --- |
| 1 | `C`, `D` |
| more than 1 | `A`, `C`, `D` |

`B` is not modified. If the same original system must be solved again, copy the
original coefficients and right-hand side back to the device arrays before the
next call. Reusing `plan` is supported as long as `nsys`, communicator layout,
and thread configuration remain compatible.

## Size and decomposition constraints

- `nsys` must be positive.
- For more than one MPI rank, every rank-local system must contain at least two
  rows because the modified local solver reads the first two rows.
- The internal reduced-system distribution requires `nsys >= nranks` so that
  every rank owns at least one reduced system.
- The helper `para` uses balanced contiguous partitioning. A global row length
  that is not divisible by `nranks` produces local lengths of either
  `floor(global_nrow / nranks)` or `ceil(global_nrow / nranks)`.

For a global row length `n3`, `n3 >= 2 * nranks` is therefore required by the
multi-rank row-length constraint.

## Algorithm flow

For one rank, `pascal_solver` launches the standard many-system Thomas solver.

For multiple ranks, it performs:

```text
modified local TDMA
  -> pack A/C/D reduced rows
  -> CUDA-aware MPI_Alltoallv
  -> assemble and solve transformed reduced systems
  -> redistribute interface solutions with MPI_Alltoallv
  -> update each full rank-local row
```

The regular solver is asynchronous except where communication requires device
completion. The profiled entry point introduces additional device
synchronization to attribute time to individual phases; it should not be
treated as having identical overhead to `pascal_solver`.

## Integrating the library

1. Build `lib/libpascal_tdma.a` and `include/pascal_tdma_cuda.mod`.
2. Compile the application with the module include path.
3. Link the static library with the same CUDA Fortran and MPI toolchain.
4. Allocate and initialize the four device arrays in the layout above.
5. Create one plan, perform compatible solves, clean the plan, and then
   finalize MPI.

The example Makefile demonstrates the required compile and link flags.

## License and citation

This component is distributed under the repository's [MIT License](../LICENSE).
Please cite the PaScaL_TDMA 2.1 paper described in the
[root citation section](../README.md#citation).
