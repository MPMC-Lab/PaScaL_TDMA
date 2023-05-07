!======================================================================================================================
!> @file        global.f90
!> @brief       This file contains a module of global parameters for the example problem of PaScaL_TDMA.
!> @details     The target example problem is the three-dimensional(3D) time-dependent heat conduction problem 
!>              in a unit cube domain applied with the boundary conditions of vertically constant temperature 
!>              and horizontally periodic boundaries.
!>
!> @author      
!>              - Ki-Ha Kim (k-kiha@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
!>              - Ji-Hoon Kang (jhkang@kisti.re.kr), Korea Institute of Science and Technology Information
!>              - Jung-Il Choi (jic@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
!>
!> @date        May 2023
!> @version     2.0
!> @par         Copyright
!>              Copyright (c) 2019-2023 Ki-Ha Kim and Jung-Il choi, Yonsei University and 
!>              Ji-Hoon Kang, Korea Institute of Science and Technology Information, All rights reserved.
!> @par         License     
!>              This project is release under the terms of the MIT License (see LICENSE file).
!======================================================================================================================

!>
!> @brief       Module for global parameters.
!> @details     This global module has simulation parameters and a subroutine to initialize the parameters. 
!>
module global

    implicit none
    double precision, parameter :: PI = acos(-1.d0)
    
    ! Physical parameters
    double precision :: Pr  !< Prandtl number
    double precision :: Ra  !< Reyleigh number

    ! Iteration steps
    integer :: Tmax             !< Maximum number of iteration steps
    
    ! Computational size for the physical domain and time discretization
    !> @{ Grid numbers in each direction
    integer :: nx,ny,nz 
    !> @}
    !> @{ Grid numbers minus 1 in each direction
    integer :: nxm,nym,nzm
    !> @}
    !> @{ Grid numbers plus 1 in each direction
    integer :: nxp,nyp,nzp
    !> @}
    double precision :: dt                  !< Length of time step
    double precision :: dtStart             !< Initial dt
    double precision :: tStart              !< Initial simulation time
    
    ! Physical size of the computational domain
    !> @{ Lengths of the physical domain
    double precision :: lx, ly, lz
    !> @}
    !> @{ Discretized grid lengths of the physical domain
    double precision :: dx, dy, dz
    !> @}
    
    ! Boundary conditions of the hot and cold walls and related parameters
    double precision :: theta_cold              !< Boundary temperature of cold wall
    double precision :: theta_hot               !< Boundary temperature of hot wall
    double precision :: alphaG                  !< Thermal expansion coefficient x gravitational acceleration
    double precision :: nu                      !< Kinematic viscosity
    double precision :: Ct                      !< Thermal diffusivity

#ifdef _CUDA
	! Thread for cuda
    integer :: thread_in_x, thread_in_y, thread_in_z, thread_in_x_pascal, thread_in_y_pascal
#endif
    
    contains 
    !>
    !> @brief       Assign global parameters.
    !> @param       np_dim      Number of MPI processes in 3D topology
    !>
    subroutine global_inputpara(np_dim)
        implicit none
        integer, intent(out) :: np_dim(0:2)

        integer :: npx, npy, npz   ! Variables to read number of processes in 3D topology
        integer :: arg_cnt
        character(len=32)   :: arg_str

        ! Namelist variables for file input
        namelist /meshes/ nx, ny, nz
        namelist /procs/ npx, npy, npz
        namelist /time/ tmax
#ifdef _CUDA
        namelist /threads/ thread_in_x, thread_in_y, thread_in_z, thread_in_x_pascal, thread_in_y_pascal
#endif

		arg_cnt = command_argument_count()

		if(arg_cnt.ne.1) then
            print *, 'Input file name is not defined. Usage:"mpirun -np number exe_file input_file" '
            stop
        endif
        call get_command_argument(arg_cnt, arg_str)

        open (unit = 1, file = arg_str)
            read (1, meshes)
            read (1, procs)
            read (1, time)
#ifdef _CUDA
			read (1, threads)
#endif
			close (1)

        np_dim(0) = npx
        np_dim(1) = npy
        np_dim(2) = npz

        ! Physical parameters
        Pr = 5.0d0; Ra = 2.d+2

        ! Computational size for the physical domain and time discretization
        nx = nx+1; ny = ny+1; nz = nz+1
        nxm = nx-1; nym = ny-1; nzm = nz-1
        nxp = nx+1; nyp = ny+1; nzp = nz+1

        dtStart = 5.0D-3; tStart = 0.d0

        ! Physical size of the computational domain
        lx = 1.0d0; ly = 1.0d0; lz = 1.0d0

        ! Boundary conditions of the hot and cold walls and related parameters
        theta_cold = -1.d0; theta_hot = 2.d0 + theta_cold
        alphaG = 1.d0; nu = 1.d0/sqrt(Ra/(alphaG*Pr*ly**3.*(theta_hot-theta_cold)))
        Ct = nu/Pr
        
    end subroutine global_inputpara

end module global