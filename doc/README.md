# PaScaL_TDMA 2.1 documentation

This directory contains the maintained manual documentation for the v2.1
source tree.

- [CUDA Fortran implementation](fortran.md)
- [CUDA C++ implementation](cuda_cxx.md)
- [CUDA C++ porting notes](cuda_cxx_porting_notes.md)
- [v2.1 release notes](release_notes/v2.1.md)
- [source and release provenance](PROVENANCE.md)
- [Study workflow and evidence boundary](study/README.md)
- [CUDA C++ engineering Study report](study/report_cuda_cxx_porting_tdma_solver_study.md)

Generated Doxygen HTML, LaTeX, and man pages from v2.0 were removed because
they documented source files and APIs that are not part of v2.1. They remain
available in the immutable [`v2.0` tag](https://github.com/MPMC-Lab/PaScaL_TDMA/tree/v2.0/doc).

`tool/gpu_power_monitor.py` is retained unchanged from v2.0. Its `@version`
header identifies that utility's original version, not the current solver
release.
