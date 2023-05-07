# PaScaL_TDMA 2.0

Parallel and Scalable Library for TriDiagonal Matrix Algorithm

PaScaL_TDMA provides an efficient and scalable computational procedure to solve many tridiagonal systems in multi-dimensional partial differential equations. The modified Thomas algorithm proposed by Laszlo et al.(2016) and the newly designed communication scheme have been used to reduce the communication overhead in solving many tridiagonal systems.

This library is for both single and many tridiagonal systems of equations. The main algorithm for a tridiagonal matrix consists of the following five steps: 

- (1) Transform the partitioned submatrices in the tridiagonal systems into modified submatrices:
        Each computing core transforms the partitioned submatrices in the tridiagonal systems of equations into the modified forms by applying modified Thomas algorithm.
- (2) Construct reduced tridiagonal systems from the modified submatrices:
        The reduced tridiagonal systems are constructed by collecting the first and last rows of the modified submatrices from each core using MPI_Ialltoallw.
- (3) Solve the reduced tridiagonal systems:
        The reduced tridiagonal systems constructed in Step 2 are solved by applying the Thomas algorithm.
- (4) Distribute the solutions of the reduced tridiagonal systems:
        The solutions of the reduced tridiagonal systems in Step 3 are distributed to each core using MPI_Ialltoallw.
        This communication is an exact inverse of the communication in Step 2.
- (5) Update the other unknowns in the modified tridiagonal systems:
        The remaining unknowns in the modified submatrices in Step 1 are solved in each computing core with the solutions obtained in Step 3 and Step 4.
    
Step 1 and Step 5 are similar to the method proposed by Laszlo et al.(2016) which uses parallel cyclic reduction (PCR) algorithm to build and solve the reduced tridiagonal systems. Instead of using the PCR, we develop an all-to-all communication scheme using the MPI_Ialltoall function after the modified Thomas algorithm is executed. The number of coefficients for the reduced tridiagonal systems are greatly reduced, so we can avoid the communication bandwidth problem, which is a main bottleneck for all-to-all communications. Our algorithm is also distinguished from the work of Mattor et al. (1995) which assembles the undetermined coefficients of the temporary solutions in a single processor using MPI_Gather, where load imbalances are serious.


# CUDA implementation in PaScaL_TDMA 2.0
In PaScaL_TDMA 2.0, multi-GPU acceleration is implemented using NVIDIA CUDA. CUDA-related features are as follows:
- (1) Incorporation of CUDA kernels into the loop structures of the existing algorithm, that are modified to exploit more GPU threads.
- (2) Utilization of shared memory using pipeline copy of variables in device memory to reduce the amount of device memory access.
- (3) CUDA-aware MPI communication for rapid communication with the support of hardward
- (4) Use of 3D-array for practical applications and accordingly the use of CUDA threads more than in a single dimension with the 3-D array.
- (5) Depreciation on functions for single tridiagonal matrix as they are rarely used for three-dimensional problems.


# Authors
- Ki-Ha Kim (k-kiha@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University (v1.0, v2.0)
- Mingyu Yang (yang926@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University (v2.0)
- Ji-Hoon Kang (jhkang@kisti.re.kr), Korea Institute of Science and Technology Information (v1.0, v2.0)
- Oh-Kyoung Kwon (okkwon@kisti.re.kr), Korea Institute of Science and Technology Information (v2.0)
- Jung-Il Choi (jic@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University (v1.0, v2.0)

# Usage
## Downloading PaScaL_TDMA
The repository can be cloned as follows:

```
git clone https://github.com/MPMC-Lab/PaScaL_TDMA.git
```
Alternatively, the source files can be downloaded through github menu 'Download ZIP'.

## Compile
### Prerequisites
Prerequisites to compile PaScaL_TDMA are as follows:
* MPI
* fortran compiler (`nvfortran` for GPU runs, NVIDIA HPC SKD 21.1 or higher)

### Compile and build
* Build PaScaL_TDMA
    ```
	make lib
	```
* Build an example problem after build PaScaL_TDMA

    ```
	make example
	```
* Build all

    ```
	make all
	```
### Mores on compile option
The `Makefile` in root directory is to compile the source code, and is expected to work for most systems. The 'Makefile.inc' file in the root directory can be used to change the compiler (and MPI wrapper) and a few pre-defined compile options depending on compiler, execution environment and et al.

## Running the example
After building the example file, an executable binary, `a.out`, is built in the `run` folder. The `PARA_INPUT.inp` file in the `run` folder is a pre-defined input file, and the `a.out` can be executed as follows:
    ```
	mpirun -np 8 ./a.out ./PARA_INPUT.inp
    ```
## GPU power monitoring
In the `tool` folder, there is a Python script `gpu_power_monitor.py` that can be used to monitor and print real-time GPU power usage with timestamps. To use this script, you will need to install the `pynvml` library.

# Folder structure
* `src` : source files of PaScaL_TDMA 2.0.
* `example` : source files of an example problem for 3D heat-transfer equation.
* `include` : header files are created after building
* `lib` : a static library of PaScaL_TDMA 2.0 is are created after building
* `doc` : documentation
* `run` : an executable binary file for the example problem is created after building.
* `tool` : contains useful scripts and tools.

# Cite
Please use the following bibtex, when you refer to this project.

    @article{kkpc2020,
        title = "PaScaL_TDMA: A library of parallel and scalable solvers for massive tridiagonal system",
        author = "Kim, Ki-Ha and Kang, Ji-Hoon and Pan, Xiaomin and Choi, Jung-Il",
        journal = "Computer Physics Communications",
        volume = "260",
        pages = "107722",
        year = "2021",
        issn = "0010-4655",
        doi = "https://doi.org/10.1016/j.cpc.2020.107722"
    }

    @misc{PaScaL_TDMA2019,
        title  = "Parallel and Scalable Library for TriDiagonal Matrix Algorithm",
        author = "Kim, Ki-Ha and Kang, Ji-Hoon and Choi, Jung-Il",
        url    = "https://github.com/MPMC-Lab/PaScaL_TDMA",
        year   = "2019"
    }


# References
For more information, please the reference paper and [School of Mathematics and Computing (Computational Science and Engineering)](https://www.mpmc.yonsei.ac.kr/)
