# CUDA C++ porting notes and implementation contract

- Initial port date: 2026-06-30
- Initial C++ implementation commit:
  `05e36c6bed3014f506868c4d6f973c40cbdd8203`
- Fortran source: `src/pascal_tdma_cuda.f90`
- C++ implementation: `src/pascal_tdma_cuda.cu`

Fortran source SHA-256 at the 2026-07-27 documentation review:
`141bae5a1945f5009614ce299fcb262e40cd29fd396e0b54633cfd005c9ec726`

This document records the design contract used for the CUDA C++ port. It is a
public implementation note, not a claim that the two language versions are
textually identical. The current public instructions are documented in the
[root README](../README.md).

## 1. Porting goals

The first port preserves the established numerical and communication structure
while translating language and runtime mechanisms:

- CUDA Fortran kernels to CUDA C++ kernels;
- Fortran device arrays to explicit CUDA allocations and flat indexing;
- the MPI Fortran interface to the MPI C interface;
- the Fortran plan type to a move-only RAII C++ plan;
- the original device-buffer MPI path, plus an optional host-staging fallback.

The following are intentionally outside this port's scope:

- replacing Thomas/modified-Thomas with PCR, CR, or another algorithm;
- changing the two-interface-row reduction;
- changing the two `MPI_Alltoallv` communication stages;
- changing precision or supporting mixed precision;
- redesigning the data layout for a C/C++ row-major application;
- adding line-internal parallelism or unrelated performance optimizations.

## 2. Preserved solver flow

For one MPI rank:

```text
tdma_many
```

For multiple MPI ranks:

```text
modified local TDMA
  -> form two reduced rows per local system
  -> pack reduced A/C/D
  -> MPI_Alltoallv
  -> unpack transformed systems
  -> initialize transformed B
  -> solve the transformed reduced systems
  -> pack their solutions
  -> MPI_Alltoallv
  -> unpack interface solutions
  -> update full rank-local rows
```

This flow is represented by `solve_impl` and the CUDA kernels in
`src/pascal_tdma_cuda.cu`.

## 3. Layout mapping

The Fortran array declaration

```fortran
real(8), device :: A(0:nsys-1,0:nrow-1)
```

has a contiguous first index. The equivalent CUDA C++ flat offset is:

```cpp
index2(sys, row, nsys) = sys + row * nsys;
```

For the z-direction example:

```text
sys = i + j * n1
row = k
```

This mapping is a numerical interoperability requirement. A caller whose native
layout differs must transform its data or supply a compatible solver view.

## 4. Plan ownership and lifecycle

`PascalTdmaPlan` owns:

- a duplicate of the supplied MPI communicator;
- rank, rank-count, and problem-size metadata;
- balanced reduced-system partition descriptors;
- reduced arrays `ard`, `brd`, `crd`, and `drd`;
- transformed arrays `atr`, `btr`, `ctr`, and `dtr`;
- all-to-all counts and displacements;
- packed device communication buffers;
- CUDA launch configurations and the selected MPI buffer mode.

The type is non-copyable and movable. `create` first destroys any existing
state, then duplicates the communicator and allocates the work buffers.
`destroy` releases CUDA storage and, while MPI is still active, frees the
duplicated communicator. Applications should therefore destroy plans before
`MPI_Finalize`.

All ranks in the communicator must create compatible plans with the same
`nsys` and participate in the same collective solve sequence.

## 5. MPI communication modes

The primary path preserves CUDA-aware MPI behavior:

```text
device buffer -> MPI_Alltoallv -> device buffer
```

The C++ port also supports:

```text
device -> host staging -> MPI_Alltoallv -> host staging -> device
```

Mode selection occurs when the plan is created:

- unset `PASCAL_TDMA_MPI_MODE`, or any value other than `host`:
  `MpiBufferMode::DeviceDirect`;
- exact value `PASCAL_TDMA_MPI_MODE=host`:
  `MpiBufferMode::HostStaging`.

The solver synchronizes the selected CUDA stream before communication so MPI
does not observe incomplete pack kernels. The profiled entry point adds more
synchronization at phase boundaries for measurement.

## 6. Deliberate language adaptations

### 6.1 Pack and unpack descriptors

The Fortran kernels receive rank descriptors from device arrays. The C++ host
loop already has the same per-rank extents and starts, so each C++ kernel launch
receives scalar `sub0`, `sub1`, `start0`, and `start1` values. The packed region
and index mapping remain the same.

### 6.2 Shared communication storage

The first exchange carries three reduced arrays (`A`, `C`, and `D`), while the
return exchange carries only the solved `D`. The C++ plan retains both the
three-array counts (`big_counts_*`) and one-array counts (`counts_*`) and
reuses the allocated buffers, matching the two exchange extents without an
extra allocation.

### 6.3 Resource and error handling

CUDA runtime calls are checked through `PASCAL_TDMA_CUDA_CHECK`, which throws
`CudaError`. Invalid plan state and detected size mismatches throw standard C++
exceptions. The examples catch these errors and abort the MPI job so that one
failed rank does not leave peers waiting in a collective.

### 6.4 Thread configuration

The plan accepts separate block sizes for the rank-local and transformed
reduced TDMA kernels. Pack dimensions are derived from the balanced system
partition. The default public values are 128 threads for both TDMA stages.

## 7. In-place numerical contract

The C++ port preserves the Fortran solver's use of caller arrays as workspace:

| Execution path | Modified caller arrays |
| --- | --- |
| one rank | `C`, `D` |
| multiple ranks | `A`, `C`, `D` |

`B` remains unchanged. A repeated solve of the same mathematical system must
restore the original coefficient and right-hand-side arrays before each call.
This is part of the public contract rather than an implementation accident.

For multiple ranks, the local modified solver requires `nrow >= 2`, and the
reduced-system partition requires `nsys >= nranks`. The plan is created for a
fixed `nsys`; changing it requires recreating the plan.

## 8. Verification strategy and current boundary

The port was checked through:

1. source-flow comparison against the CUDA Fortran implementation;
2. identical system-contiguous layout in both examples;
3. single-rank and multi-rank execution;
4. expected-solution checks for both device-direct and host-staging modes;
5. phase-level timing and rank-wise maximum reductions;
6. multi-GPU runs on an NVIDIA H200 system.

The current curated Study dataset contains 25 CUDA C++ correctness cases. The
maximum iteration-0 absolute error is `8.860e-12`. It does not yet contain
Fortran timing rows.

The profiling drivers initialize inputs once and repeat the in-place solver.
Therefore, their iterations after iteration 0 exercise successively modified
arrays. Current iterations 1–9 timing can be used to inspect that execution
protocol and its scaling, but it is not timing for repeated independent solves
of the same original system. See the [Study documentation](study/README.md).

## 9. Maintenance rule

Changes to the port should preserve, or explicitly document a change to, all of
the following:

- solver stage ordering;
- flat array layout;
- collective count and displacement semantics;
- caller-array mutation behavior;
- plan ownership and MPI lifecycle;
- device-direct and host-staging mode selection;
- validation scope and the exact protocol behind reported measurements.
