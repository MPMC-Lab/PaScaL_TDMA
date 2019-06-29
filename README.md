# PaScaL_TDMA

Parallel and Scalable Library for Tri-Diagonal Matrix Algorithm

PaScal_TDMA solves a single or many tridiagonal systems of equations in a parallel manner using a noble all-to-all communication scheme with the modified Thomas algorithm by Laszlo et al(2016).

This library is for both of a single and many tridiagonal systems of equations. The main algorithm for a tridiagonal matrix proposed in this library consists of the following five steps: 

- (1) Transform the original partitioned tridiagonal systems of equations 
        The original partitioned tridiagonal system of equations is transformed to the modified tridiagonal system using the modified Thomas algorithm.
- (2) Build the reduced tridiagonal systems of equations
        The first and last rows of the modified tridiagonal systems are assembled by executing the proposed all-to-all communication scheme and the reduced tridiagonal systems are built.
- (3) Solve the reduced tridiagonal systems
        The reduced tridiagonal systems are solved by applying Thomas algorithm.
- (4) Distribute the solutions of the reduced tridiagonal systems
        The solutions in Step 3 are distributed to each computing cores by executing the inverse of the all-to-all communication in Step 2.
- (5) Update the solutions in the modified tridiagonal systems
        The remaining unknows are updated using the modified tridiagonal systems with the solutions obtained in Step 3 and Step 4.
    
Step 1 and Step 5 are similar to the method proposed by Laszlo, Gilles and Appleyard(2016)
which used the parallel cyclic reduction (PCR) to build and solve the reduced tridiagonal systems.
Instead of using the PCR, we develop an all-to-all communication scheme using a MPI_Ialltoall
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
For more information, please the reference paper [Paper name] and [Multi-Physics Modeling and Computation Lab.](https://www.mpmc.yonsei.ac.kr/)