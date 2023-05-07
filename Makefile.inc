###############################################################################
# Compile setting with intel compiler and intel MPI.
# Modify flags according to the system and environment
###############################################################################

# FC = mpiifort

# FLAG  = -O3 -fPIC -cpp
# FLAG += -g
# FLAG += -fcheck=all -warn-all 
# FLAG += -fp-model=precise 
# FLAG += -check bounds
# FLAG += -xMIC-AVX512

###############################################################################
# Compile setting with nvfortran and non-intel MPI.
# Modify flags according to the system and environment
###############################################################################

FC = mpifort

FLAG  = -O3 -fPIC -cpp 
FLAG += -tp=skylake 

###############################################################################
# To use CUDA
# PaScaL_TDMA is compiled for CUDA also.
###############################################################################

USE_CUDA = yes
FLAG += -mp:multicore -Mcuda -Minfo=accel 
#FLAG += -DUSE_NVTX


