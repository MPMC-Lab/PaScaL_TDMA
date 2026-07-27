# PaScaL_TDMA 2.1

Parallel and Scalable Library for TriDiagonal Matrix Algorithm

PaScaL_TDMA solves large batches of tridiagonal systems whose row direction
is distributed across MPI ranks. Version 2.1 supersedes version 2.0 as the
GPU-oriented release while preserving the five-step numerical structure of
the original PaScaL_TDMA method.

The PaScaL_TDMA 2.1 publication and program archive describe the CUDA Fortran
implementation. This repository also provides a CUDA C++17 port and
profiling/study tools as post-publication repository extensions.

## Version lineage

| Release | Main scope | Publication |
| --- | --- | --- |
| `v1.0` | Original CPU implementation | [Computer Physics Communications 260 (2021) 107722](https://doi.org/10.1016/j.cpc.2020.107722) |
| `v2.0` | Multi-GPU CUDA Fortran implementation with shared-memory buffering | [Computer Physics Communications 290 (2023) 108785](https://doi.org/10.1016/j.cpc.2023.108785) |
| `v2.1` | Register-resident CUDA Fortran implementation with a redesigned kernel interface and consolidated communication | [Computer Physics Communications 323 (2026) 110120](https://doi.org/10.1016/j.cpc.2026.110120) |

See the [changelog](CHANGELOG.md), [v2.1 release notes](doc/release_notes/v2.1.md),
and [source provenance](doc/PROVENANCE.md) for the detailed boundary between
the published implementation and later repository extensions.

## Algorithm

For a distributed solve, each rank owns a contiguous segment of every
tridiagonal line. PaScaL_TDMA performs:

1. a modified Thomas reduction on each rank-local line segment;
2. packing and exchange of the first and last reduced rows;
3. solution of the assembled reduced tridiagonal systems;
4. redistribution of the interface solutions; and
5. reconstruction of each full rank-local solution.

The current source uses two collective phases per distributed solve. The
forward phase packs reduced `A`, `C`, and `D` fields into one
`MPI_Alltoallv` exchange. After the reduced solve, the backward phase returns
interface `D` values through a second `MPI_Alltoallv`. The publication's
description of a consolidated all-to-all refers to combining the coefficient
fields used for reduced-system assembly; it does not remove the later
interface-solution exchange.

For one MPI rank, the implementation directly applies the many-system Thomas
solver.

## What changed in version 2.1

Compared with PaScaL_TDMA 2.0, version 2.1:

- replaces shared-memory buffering with register-resident elimination
  kernels;
- maps adjacent CUDA threads to adjacent tridiagonal lines for coalesced
  global-memory access;
- replaces assumed-shape kernel arrays with fixed-bound arrays and explicit
  scalar arguments, avoiding implicit descriptor transfers at launch; and
- consolidates the fields used for reduced-system assembly into a packed
  all-to-all exchange while reconstructing the unit diagonal locally.

These changes preserve the numerical reduction and reconstruction sequence.
Performance results in the article are specific to the hardware, problem
sizes, and placement used in that study.

## Implementations

| Implementation | Source | Scope |
| --- | --- | --- |
| CUDA Fortran | `src/pascal_tdma_cuda.f90` | Published PaScaL_TDMA 2.1 implementation lineage; CUDA-aware MPI device buffers |
| CUDA C++17 | `src/pascal_tdma_cuda.cu`, `src/pascal_tdma_cuda.hpp` | Repository port of the same solver flow; device-direct MPI by default and optional host staging |

Both implementations solve non-cyclic, double-precision systems. The CUDA C++
implementation was not the programming-language implementation evaluated in
the PaScaL_TDMA 2.1 article.

## Requirements

- NVIDIA GPU and compatible driver;
- NVIDIA HPC SDK with CUDA Fortran support;
- CUDA Toolkit with `nvcc` for the CUDA C++ port;
- MPI Fortran and C++ compiler wrappers;
- CUDA-aware MPI for CUDA Fortran and for the default CUDA C++ device-buffer
  path; and
- GNU Make.

`FC` must select an MPI wrapper backed by NVIDIA `nvfortran`. The checked-in
default is `FC=mpifort`; override the wrapper and `CUDA_ARCH` for the target
system. Only the CUDA C++ port supplies an explicit host-staging fallback.

## Build

The historical MPMC entry points retain their scope:

```bash
make lib CUDA_ARCH=90 FC=mpifort
make example CUDA_ARCH=90 FC=mpifort
make all CUDA_ARCH=90 FC=mpifort
```

| Target | Output |
| --- | --- |
| `make lib` | `lib/libpascal_tdma.a` and Fortran modules in `include/` |
| `make example` | the Fortran library and `run/a.out` |
| `make all` | the same library and Fortran example |
| `make cuda-cxx` | `lib/libpascal_tdma_cuda.a`, the public header in `include/`, and two C++ examples |
| `make cuda-cxx-profile` | the C++ library and `run/ex_tdma_profile` |
| `make study` | both libraries and both corresponding profile drivers |
| `make all-implementations` | every library, example, and profile driver |

For the added implementations:

```bash
make cuda-cxx CUDA_ARCH=90 MPICXX=mpicxx
make study CUDA_ARCH=90 FC=mpifort MPICXX=mpicxx
```

Compiler and architecture defaults are centralized in `Makefile.inc`.
`make clean` removes generated build products while retaining `run/job.sh` and
the checked-in Study data.

## Run

CUDA Fortran example:

```bash
mpirun -np 4 ./run/a.out
```

CUDA C++ example with CUDA-aware MPI:

```bash
mpirun -np 4 ./run/ex_tdma_zdirection
```

CUDA C++ host-staging fallback:

```bash
PASCAL_TDMA_MPI_MODE=host \
  mpirun -np 4 ./run/ex_tdma_zdirection
```

Only the exact value `host` selects host staging. The examples assign a GPU
using `rank % visible_device_count`; configure rank placement to avoid
unintended GPU oversubscription.

## Solver data contract

Both implementations use system-contiguous storage:

```text
offset(sys, row) = sys + row * nsys
```

`D` contains the right-hand side on entry and the solution on return. A
one-rank solve modifies `C` and `D`; a distributed solve modifies `A`, `C`, and
`D`. `B` is not modified. Restore the original modified arrays before solving
the same mathematical system again.

For multiple ranks, each rank must have `nrow >= 2`, and the current reduced
system partition requires `nsys >= number of ranks`.

## Repository layout

- `src/`: CUDA Fortran and CUDA C++ solver sources;
- `examples/`: standalone examples and corresponding profile drivers;
- `run/`: generated executables and the cluster job example;
- `include/`, `lib/`: generated headers/modules and static libraries;
- `doc/`: manual documentation, release notes, provenance, and Study data;
- `tool/`: profiling, environment, cleanup, and GPU power utilities.

The [documentation index](doc/README.md) links the implementation manuals and
Study report. Generated v2.0 Doxygen pages were removed because they described
source and APIs no longer present; the historical pages remain available from
the immutable `v2.0` tag.

## Study boundary

Build and preview the repository Study workflow with:

```bash
make study CUDA_ARCH=90 FC=mpifort MPICXX=mpicxx
DRY_RUN=1 STUDY_PRESET=quick ./tool/run_study_sweep.sh
```

The current curated dataset contains 25 CUDA C++ cases and no CUDA Fortran
timing rows. Correctness is recorded after iteration 0. Later timed iterations
reuse arrays transformed in place by earlier solver calls, so the data do not
establish repeated independent-solve timing or a Fortran-versus-C++ comparison.
See the [Study documentation](doc/study/README.md) before interpreting the
tables.

## Release contributors

The version labels below record public release attribution; per-file
copyright and article authorship are documented separately.

- Ki-Ha Kim (`v1.0`, `v2.0`, `v2.1`)
- Mingyu Yang (`v2.0`)
- Ji-Hoon Kang (`v1.0`, `v2.0`, `v2.1`)
- Oh-Kyoung Kwon (`v2.0`)
- Jung-Il Choi (`v1.0`, `v2.0`, `v2.1`)
- Dongjin Lee (`v2.1`)
- Junhwan Lee (`v2.1`)
- Sehyeong Oh (`v2.1`)
- Seungwon Lee (`v2.1`)

## Citation

Use [CITATION.cff](CITATION.cff) for machine-readable metadata. Cite the
PaScaL_TDMA 2.1 article for the published CUDA Fortran release and the earlier
articles when discussing the original or version 2.0 methods. The publication
snapshot is archived as
[Mendeley Data version 3](https://doi.org/10.17632/49z6fh94z3.3); the current
repository also contains post-archive fixes, profiling support, and a CUDA C++
port.

## License

PaScaL_TDMA is distributed under the [MIT License](LICENSE). Historical and
imported-code copyright notices are preserved in that file.
