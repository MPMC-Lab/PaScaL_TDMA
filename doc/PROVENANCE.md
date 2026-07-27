# Source and release provenance

This document records the repository lineage and the boundary between the
PaScaL_TDMA 2.1 publication snapshot and the maintained source tree. It does
not replace the release notes or the software license.

## MPMC release lineage

The v2.1 work is based on the existing `MPMC-Lab/PaScaL_TDMA` history. The
previous release tags are preserved and are not recreated or rewritten.

| Reference | Commit | Date | Note |
| --- | --- | --- | --- |
| `v1.0` | `f202864854909c40a4254372d4223b77eada86b4` | 2021-06-14 | Original CPU release line |
| `v2.0` | `b5d15c72cb4386b6a0bf29eada7b396bef630024` | 2023-05-07 | CUDA-enabled v2.0 release line |
| pre-v2.1 `master` | `e68f2b857bf730faf609eecad7ed710c41fe64e7` | 2023-05-14 | v2.0 plus a README update |

The v2.1 integration branch was created from the pre-v2.1 `master`, so a v2.1
merge and tag continue the existing v1.0 -> v2.0 -> v2.1 Git history.

## PaScaL_TDMA 2.1 publication snapshot

The v2.1 method is described in:

> Ki-Ha Kim, Dongjin Lee, Junhwan Lee, Sehyeong Oh, Seungwon Lee, Ji-Hoon
> Kang, and Jung-Il Choi, “PaScaL_TDMA 2.1: A register-resident multi-GPU
> tridiagonal matrix solver with optimized communication for large-scale CFD
> simulations,” *Computer Physics Communications* 323 (2026) 110120.
> <https://doi.org/10.1016/j.cpc.2026.110120>

The publication snapshot is archived as Mendeley Data version 3:
<https://doi.org/10.17632/49z6fh94z3.3>.

This Git repository is the maintained development version. It is not claimed
to be byte-identical to the archived snapshot.

## Integration source

The maintained v2.1 implementation was imported from
[`k-kiha/PaScaL_TDMAcuda`](https://github.com/k-kiha/PaScaL_TDMAcuda) at commit
[`e637d5ecab0f08308eea83fbbb1872ead7ae07c5`](https://github.com/k-kiha/PaScaL_TDMAcuda/commit/e637d5ecab0f08308eea83fbbb1872ead7ae07c5).

Only file contents were imported. The independent repository's Git history
was not grafted into the MPMC history; the external commit links below retain
that provenance.

| MPMC v2.1 path | Integration-source path | SHA-256 |
| --- | --- | --- |
| `src/pascal_tdma_cuda.f90` | `Fortran_Original/src/PaScaL_TDMA_cuda.f90` | `141bae5a1945f5009614ce299fcb262e40cd29fd396e0b54633cfd005c9ec726` |
| `src/pascal_tdma_cuda.cu` | `CUDA_CXX_Port/src/pascal_tdma_cuda.cu` | `0d5b0543dcefe0fa101a5a4b16cb3a9d4f44743151182c9a73af4250dff4e6c3` |
| `src/pascal_tdma_cuda.hpp` | `CUDA_CXX_Port/include/pascal_tdma_cuda.hpp` | `17a5dc5ac2ecc070dcffcfe378a42a49b78f310282fcf1aeaad64c648bc23fc1` |
| `examples/ex_tdma_zdirection.f90` | `Fortran_Original/examples/ex_tdma_zdirection.f90` | `064d1d45a696196beedeeab292af676f5eb95fe0dcc3e07933c60afb0e4db636` |
| `examples/ex_tdma_zdirection.cu` | `CUDA_CXX_Port/examples/ex_tdma_zdirection.cu` | `40a632bc1b61fe3343dcd68b59cf71daa04ea6434b9f3e6eb50298744be8ba35` |
| `examples/ex_tdma_profile.cu` | `CUDA_CXX_Port/examples/ex_tdma_profile.cu` | `bc87e1a06293f1ba6e1d60d4b84743f9035706f3ddbcc4145f8790d2db6b861d` |
| `examples/example_fortran_profile.f90` | `Study/example_fortran_profile.f90` | `92061affac3c34d453069ab678578170126c898313682008cc88c870fe8653c1` |
| `examples/example_cuda_cxx_profile.cu` | `Study/example_cuda_cxx_profile.cu` | `09125c50b18250536a796cc70742cc57af1d44352d71ed88ae108eb18c1390de` |

The path layout was adapted to the MPMC repository structure. The hashes above
record that the numerical solver and example/profile source files were
unchanged during that path adaptation.

## Maintained changes beyond the publication snapshot

The integration source contains maintenance and repository work beyond the
archived publication snapshot, including:

- a CUDA Fortran request-handling fix in
  [`d69ae952`](https://github.com/k-kiha/PaScaL_TDMAcuda/commit/d69ae95209b69752504ebe3888721d759ae3020b);
- the CUDA C++ implementation introduced in
  [`05e36c6b`](https://github.com/k-kiha/PaScaL_TDMAcuda/commit/05e36c6bed3014f506868c4d6f973c40cbdd8203);
- profiling and Study support developed through
  [`902705ea`](https://github.com/k-kiha/PaScaL_TDMAcuda/commit/902705ea45ddb3d6fe6e20de7ee976a957e50a24)
  and later commits.

The paper and publication archive describe the CUDA Fortran implementation.
The CUDA C++ port and the repository Study workflow are maintained repository
extensions and should not be presented as implementations evaluated in the
paper.

## Communication terminology

The paper describes a consolidated `MPI_Alltoall` exchange at the algorithmic
level. The maintained implementations use `MPI_Alltoallv` because message
counts can differ among ranks. A distributed solve has two collective phases:

1. packed reduced-system coefficients are exchanged for the global reduced
   solve; and
2. solved interface values are exchanged for local back substitution.

Accordingly, a consolidated all-to-all refers to consolidating coefficient
fields in the reduced-system assembly phase. It does not mean that the entire
solver performs only one collective call.

## Validation boundary

Performance values in the v2.1 article are tied to the hardware, problem
sizes, and configurations reported there. They are not general performance
guarantees for every build or machine.

The checked-in Study data are a separate repository study. At integration, the
curated dataset contains CUDA C++ measurements and no CUDA Fortran timing rows.
Correctness is recorded after the first iteration only, while later iterations
reuse arrays transformed in place by the solver. These data therefore do not
establish an independent Fortran-versus-CUDA-C++ comparison and should be
interpreted with the accompanying Study documentation.

## Release integrity

The historical `v1.0` and `v2.0` references remain immutable. The `v2.1` tag
should be created from the commit merged into the MPMC default branch, and the
GitHub release should cite that tag. Merging a pull request transfers the
proposed changes into upstream history, but it does not transfer or create a
fork-side GitHub Release.
