# PaScaL_TDMA

Parallel and Scalable Library for Tri-Diagonal Matrix Algorithm

PaScal_TDMA provides an efficient and scalable computational procedure to solve many tridiagonal systems in multi-dimensional partial differential equations. The modified Thomas algorithm by Laszlo et al.(2016) and newly designed communication scheme are used to reduce the communication overhead in solving many tridiagonal systems.

This library is for both of a single and many tridiagonal systems of equations. The main algorithm for a tridiagonal matrix proposed in this library consists of the following five steps: 

- (1) Transform the partitioned sub-matrices in the tridiagonal systems into modified sub-matrices
        Each computing core transforms the partitioned sub-matrices in the tridiagonal systems of equations into the modified forms by applying the modified Thomas algorithm.
- (2) Construct reduced tridiagonal systems from the modified sub-matrices
        The reduced tridiagonal systems are constructed by collecting the first and last row of the modified sub-matrices from each cores using MPI_Ialltoallw.
- (3) Solve the reduced tridiagonal systems
        The reduced tridiagonal systems constructed in Step 2 are solved by applying the Thomas algorithm.
- (4) Distribute the solution of reduced tridiagonal system
        The solutions of reduced tridiagonal systems in Step 3 are distributed to each core using MPI_Ialltoallw.
        This communication is an exact inverse of communication in Step 2.
- (5) Update the other unknowns in the modified tridiagonal systems
        The remaining unknowns of the modified sub-matrices in Step 1 are solved in each computing core with the solutions obtained in Step 3 and Step 4.
    
Step 1 and Step 5 are similar to the method proposed by Laszlo, Gilles and Appleyard(2016)
which used the parallel cyclic reduction (PCR) to build and solve the reduced tridiagonal systems.
Instead of using the PCR, we develop an communication scheme using a MPI_Ialltoall
function after the modified Thomas algorithm is execued. The number of coefficients for
the reduced tridiagonal systems are greatly reduced, so we can avoid the communication 
bandwidth problem which is a main bottle-neck of all-to-all communications.
Our algorithm is also distinguished from the work of Mattor, Williams, Hewette (1995) which
assembles the undetermined coefficients of the temporary solutions in a single processor 
using MPI_Gather, where load imbalances are serious.

# Authors
- Kiha Kim (k-kiha@yonsei.ac.kr), Multi-Physics Modeling and Computation Lab., Yonsei University
- Ji-Hoon Kang (jhkang@kisti.re.kr), Korea Institute of Science and Technology Information
- Jung-Il Choi (jic@yonsei.ac.kr), Multi-Physics Modeling and Computation Lab., Yonsei University

# References
For more information, please the reference paper (in preparation) and [Multi-Physics Modeling and Computation Lab.](https://www.mpmc.yonsei.ac.kr/)