# PaScaL_TDMA 2.1 top-level build
#
# The historical MPMC entry points keep their original scope:
#   make lib      CUDA Fortran library
#   make example  CUDA Fortran example (after the library)
#   make all      library and CUDA Fortran example

include Makefile.inc

.PHONY: all lib example fortran fortran-lib fortran-example \
	libs cuda-cxx cuda-cxx-lib cuda-cxx-example cuda-cxx-profile \
	study all-implementations dirs clean veryclean help

all: lib example

dirs:
	mkdir -p include lib run

lib: dirs
	$(MAKE) -C src fortran-lib

example: lib
	$(MAKE) -C examples fortran-example

fortran-lib: lib

fortran-example: example

fortran: example

libs: lib cuda-cxx-lib

cuda-cxx-lib: dirs
	$(MAKE) -C src cuda-cxx-lib

cuda-cxx-example: cuda-cxx-lib
	$(MAKE) -C examples cuda-cxx-example

cuda-cxx-profile: cuda-cxx-lib
	$(MAKE) -C examples cuda-cxx-profile

cuda-cxx: cuda-cxx-example cuda-cxx-profile

study: lib cuda-cxx-lib
	$(MAKE) -C examples study

all-implementations: all cuda-cxx study

clean:
	$(MAKE) -C examples clean
	$(MAKE) -C src clean

veryclean: clean
	@rmdir include lib 2>/dev/null || true

help:
	@echo "PaScaL_TDMA 2.1 build targets"
	@echo "  make lib                  CUDA Fortran library"
	@echo "  make example              CUDA Fortran example (run/a.out)"
	@echo "  make all                  legacy scope: lib + example"
	@echo "  make cuda-cxx             CUDA C++ library and examples"
	@echo "  make cuda-cxx-profile     CUDA C++ profiling example"
	@echo "  make study                both corresponding profile drivers"
	@echo "  make all-implementations  all libraries, examples, and drivers"
	@echo "  make clean                remove generated build artifacts"
