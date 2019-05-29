# STDMA

Scalable tri-diagonal matrix(TDM) solver for many tridiagonal systems of equations using a noble all-to-all communication scheme for the reduced unknowns based onthe modified Thomas algorithm by Laszlo, Gilles and Appleyard(2016)

This library is for both of a single and many tridiagonal systems of equations. The main algorithm for a tridiagonal matrix proposed in this library consists of the following three steps: 
    (1) Lower/Upper element reduction step : creates reduced tridiagonals system of equations.
    (2) Reduced TDMA step : Solves the reduced tridiagonal systems of equations.
    (3) Update step : Updates solutions of the origianl tridiagonal systems of equations
                      using the solutions of the reduced tridiagonal systems of equations.

Our algorithm is close to the algorithm by Mattor, Williams, Hewette (1995) except that we use MPI_Alltoall communication for reduced tridiagonal systems to distribute the many tridiagonal systems of equations almost equally to processes, while Mattor et al used MPI_Gather to gather and solve the reduced equations on the master process only.
The number of coefficient is greatly reduced by local reduction in Step(1), so we can avoid the communication bandwidth problem which happens to be common in MPI_Alltoall communication.
Our algorithm is also distinguished from the work of Laszlo, Gilles and Appleyard(2016) which used parallel cyclic reduction to solve the reduced tridiagonal systems of equations.

# Authors      
Kiha Kim (k-kiha@yonsei.ac.kr), Multi-Physics Modeling and Computation Lab., Yonsei University
Ji-Hoon Kang (jhkang@kisti.re.kr), Korea Institute of Science and Technology Information
Jung-Il Choi (jic@yonsei.ac.kr), Multi-Physics Modeling and Computation Lab., Yonsei University

# References
For more information, please the reference paper [Paper name] and [Multi-Physics Modeling and Computation Lab.](https://www.mpmc.yonsei.ac.kr/)