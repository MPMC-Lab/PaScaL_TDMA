########################################################################
# Common build configuration for PaScaL_TDMA 2.1
#
# FC must be an MPI wrapper backed by NVIDIA nvfortran. Command-line
# assignments override these portable defaults, for example:
#
#   make all FC=ftn CUDA_ARCH=90
########################################################################

FC       = mpifort
NVCC     = nvcc
MPICXX   = mpicxx
AR       = ar
RANLIB   = ranlib

CUDA_ARCH = 90
OPT       = -O3

# CUDA Fortran compile and link flags.
FLAG      = $(OPT) -fPIC -Mfree -cpp
CUDAFLAG  = -cuda -gpu=cc$(CUDA_ARCH)

# CUDA C++ compile and link flags. The MPI wrapper is used as nvcc's host
# compiler so MPI headers and libraries are supplied consistently.
GENCODE   = -gencode arch=compute_$(CUDA_ARCH),code=sm_$(CUDA_ARCH)
NVCCFLAGS = -std=c++17 $(OPT) $(GENCODE) -ccbin $(MPICXX)
