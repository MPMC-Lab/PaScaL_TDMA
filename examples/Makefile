include ../Makefile.inc

BINDIR       := ../run
INCDIR       := ../include
LIBDIR       := ../lib
FORTRAN_LIB  := $(LIBDIR)/libpascal_tdma.a
CXX_LIB      := $(LIBDIR)/libpascal_tdma_cuda.a
FORTRAN_SOURCE := ../src/pascal_tdma_cuda.f90
CXX_SOURCE     := ../src/pascal_tdma_cuda.cu
CXX_HEADER     := ../src/pascal_tdma_cuda.hpp

FORTRAN_EXE  := $(BINDIR)/a.out
CXX_EXE      := $(BINDIR)/ex_tdma_zdirection
CXX_PROFILE  := $(BINDIR)/ex_tdma_profile
STUDY_F_EXE  := $(BINDIR)/example_fortran_profile
STUDY_CXX_EXE := $(BINDIR)/example_cuda_cxx_profile

.PHONY: all fortran-example cuda-cxx-example cuda-cxx-profile study clean

all: fortran-example

fortran-example: $(FORTRAN_EXE)

cuda-cxx-example: $(CXX_EXE)

cuda-cxx-profile: $(CXX_PROFILE)

study: $(STUDY_F_EXE) $(STUDY_CXX_EXE)

$(FORTRAN_EXE): ex_tdma_zdirection.f90 $(FORTRAN_LIB) | $(BINDIR)
	$(FC) $(FLAG) $(CUDAFLAG) -I$(INCDIR) $< \
		-L$(LIBDIR) -lpascal_tdma -o $@

$(CXX_EXE): ex_tdma_zdirection.cu $(CXX_LIB) | $(BINDIR)
	$(NVCC) $(NVCCFLAGS) -I../src $< $(CXX_LIB) -o $@

$(CXX_PROFILE): ex_tdma_profile.cu $(CXX_LIB) | $(BINDIR)
	$(NVCC) $(NVCCFLAGS) -I../src $< $(CXX_LIB) -o $@

$(STUDY_F_EXE): example_fortran_profile.f90 $(FORTRAN_LIB) | $(BINDIR)
	$(FC) $(FLAG) $(CUDAFLAG) -I$(INCDIR) $< \
		-L$(LIBDIR) -lpascal_tdma -o $@

$(STUDY_CXX_EXE): example_cuda_cxx_profile.cu $(CXX_LIB) | $(BINDIR)
	$(NVCC) $(NVCCFLAGS) -I../src $< $(CXX_LIB) -o $@

$(FORTRAN_LIB): $(FORTRAN_SOURCE)
	$(MAKE) -C ../src fortran-lib

$(CXX_LIB): $(CXX_SOURCE) $(CXX_HEADER)
	$(MAKE) -C ../src cuda-cxx-lib

$(BINDIR):
	mkdir -p $@

clean:
	rm -f $(FORTRAN_EXE) $(CXX_EXE) $(CXX_PROFILE)
	rm -f $(STUDY_F_EXE) $(STUDY_CXX_EXE)
	rm -f *.o *.mod *__genmod.f90
