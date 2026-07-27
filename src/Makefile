include ../Makefile.inc

LIBDIR        := ../lib
INCDIR        := ../include
OBJDIR        := obj
FORTRAN_OBJDIR := $(OBJDIR)/fortran
CXX_OBJDIR     := $(OBJDIR)/cuda-cxx

FORTRAN_SRC   := pascal_tdma_cuda.f90
FORTRAN_OBJ   := $(FORTRAN_OBJDIR)/pascal_tdma_cuda.o
FORTRAN_LIB   := $(LIBDIR)/libpascal_tdma.a

CXX_SRC       := pascal_tdma_cuda.cu
CXX_PUBLIC_H  := pascal_tdma_cuda.hpp
CXX_OBJ       := $(CXX_OBJDIR)/pascal_tdma_cuda.o
CXX_LIB       := $(LIBDIR)/libpascal_tdma_cuda.a
CXX_STAGED_H  := $(INCDIR)/pascal_tdma_cuda.hpp

.PHONY: all fortran-lib cuda-cxx-lib clean

all: fortran-lib

fortran-lib: $(FORTRAN_LIB) | $(INCDIR)
	cp $(FORTRAN_OBJDIR)/*.mod $(INCDIR)/

cuda-cxx-lib: $(CXX_LIB) $(CXX_STAGED_H)

$(FORTRAN_LIB): $(FORTRAN_OBJ) | $(LIBDIR)
	rm -f $@
	$(AR) rcs $@ $^
	$(RANLIB) $@

$(FORTRAN_OBJ): $(FORTRAN_SRC) | $(FORTRAN_OBJDIR)
	$(FC) $(FLAG) $(CUDAFLAG) -c $< -o $@ -module $(FORTRAN_OBJDIR)

$(CXX_LIB): $(CXX_OBJ) | $(LIBDIR)
	rm -f $@
	$(AR) rcs $@ $^
	$(RANLIB) $@

$(CXX_OBJ): $(CXX_SRC) $(CXX_PUBLIC_H) | $(CXX_OBJDIR)
	$(NVCC) $(NVCCFLAGS) -I. -c $< -o $@

$(CXX_STAGED_H): $(CXX_PUBLIC_H) | $(INCDIR)
	cp $< $@

$(LIBDIR) $(INCDIR) $(FORTRAN_OBJDIR) $(CXX_OBJDIR):
	mkdir -p $@

clean:
	rm -rf $(OBJDIR)
	rm -f $(FORTRAN_LIB) $(CXX_LIB)
	rm -f $(INCDIR)/*.mod $(CXX_STAGED_H)
