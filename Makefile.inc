# Main compiler with MPI. Use mpiifort for intel MPI.
FC = mpifort

# Complie flag
# For debugging
FLAG = -O0 -g -fcheck=all --warn-all
# FLAG = -O0 -g -fPIC -check all -warn all -WB
# For moderate debugging
#FLAG = -O1 -fPIC -fp-model=precise -check bounds #-xMIC-AVX512
# For optimized run in Intel Xeon Phi
#FLAG = -O3 -fPIC -fp-model=precise -xMIC-AVX512

