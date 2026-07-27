# Changelog

All notable changes to PaScaL_TDMA are documented in this file. Historical
release tags remain unchanged.

## [2.1] - Unreleased

PaScaL_TDMA 2.1 supersedes version 2.0 as the current GPU-oriented release
while preserving the original five-step distributed TDMA method.

### Added

- Register-resident CUDA Fortran kernels for the non-cyclic many-system solver.
- A CUDA C++17 port that preserves the system-contiguous layout and distributed
  reduction/reconstruction flow.
- Optional host-staging communication for the CUDA C++ port through
  `PASCAL_TDMA_MPI_MODE=host`.
- Phase-profiled Fortran and C++ solve entry points, example drivers, sweep
  tooling, and a checked-in engineering Study.
- Machine-readable citation metadata, release notes, and explicit source
  provenance documentation.

### Changed

- Replaced the shared-memory buffering used in version 2.0 with thread-local,
  register-resident elimination kernels.
- Mapped adjacent CUDA threads to adjacent tridiagonal lines for coalesced
  global-memory access.
- Replaced the version 2.0 CUDA Fortran public API (`ptdma_plan_many_cuda` and
  `PaScaL_TDMA_*_cuda` procedures) with `ptdma_plan_cuda`,
  `pascal_plan_create`, `pascal_solver`, and `pascal_plan_clean`.
- Replaced assumed-shape kernel arrays with fixed-bound arrays and explicit
  scalar arguments in the CUDA Fortran implementation.
- Consolidated reduced `A`, `C`, and `D` fields into one packed forward
  `MPI_Alltoallv` exchange. A second `MPI_Alltoallv` returns the reduced
  interface solution before local reconstruction.
- Reorganized the release under the existing MPMC repository shell (`src/`,
  `examples/`, `run/`, `doc/`, and `tool/`).

### Removed

- The version 2.0 CPU implementation and its public API.
- The version 2.0 cyclic CUDA solve entry point; version 2.1 in this repository
  is non-cyclic.
- The version 2.0 heat-transfer example.
- Generated version 2.0 Doxygen HTML, LaTeX, and man pages, which described
  source and APIs no longer present in version 2.1.

### Compatibility and requirements

- Solver arrays use double precision and system-contiguous storage,
  `offset(sys, row) = sys + row * nsys`.
- A one-rank solve modifies `C` and `D`; a distributed solve modifies `A`, `C`,
  and `D`. `B` is not modified.
- Distributed solves require at least two local rows on every rank and
  currently require `nsys >= number of ranks`.
- CUDA Fortran requires NVIDIA GPUs and CUDA-aware MPI. The CUDA C++ port is
  also NVIDIA CUDA based, but its explicit host-staging mode does not require
  CUDA-aware MPI.

### Evidence boundary

The PaScaL_TDMA 2.1 article and Mendeley Data program archive describe a CUDA
Fortran implementation. The CUDA C++ port and the checked-in Study workflow
are later repository extensions. Current curated Study data contain CUDA C++
rows only, record correctness after iteration 0, and repeat the in-place solver
without restoring its inputs; they are not a Fortran-versus-C++ performance
comparison.

Hardware- and problem-specific performance results are reported in the
[PaScaL_TDMA 2.1 article](https://doi.org/10.1016/j.cpc.2026.110120) and are not
generalized here.

## [2.0] - 2023-05-07

- Added the first multi-GPU CUDA Fortran release.
- Used shared-memory buffering to reduce global-memory traffic.
- Added CUDA-aware MPI communication, three-dimensional examples, and
  GPU-oriented tooling.
- Published as [Computer Physics Communications 290 (2023) 108785](https://doi.org/10.1016/j.cpc.2023.108785).

## [1.0] - 2021-06-14

- Released the original CPU implementation of the parallel and scalable
  distributed TDMA method.
- Published as [Computer Physics Communications 260 (2021) 107722](https://doi.org/10.1016/j.cpc.2020.107722).

[2.0]: https://github.com/MPMC-Lab/PaScaL_TDMA/releases/tag/v2.0
[1.0]: https://github.com/MPMC-Lab/PaScaL_TDMA/releases/tag/v1.0
