!======================================================================================================================
!> @file        mpi_topology.f90
!> @brief       This file contains a module of communication topology for the example problem of PaScaL_TDMA.
!> @details     The target example problem is the three-dimensional time-dependent heat conduction problem 
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
!> @brief       Module for creating the cartesian topology of the MPI processes and subcommunicators.
!> @details     This module has three subcommunicators in each-direction and related subroutines.
!>
module mpi_topology
    use mpi
    implicit none

    integer, public :: mpi_world_cart       !< Communicator for cartesian topology
    integer, public :: np_dim(0:2)          !< Number of MPI processes in 3D topology
    logical, public :: period(0:2)          !< Periodicity in each direction

    !> @brief   Type variable for the information of 1D communicator
    type, public :: cart_comm_1d
        integer :: myrank                   !< Rank ID in current communicator
        integer :: nprocs                   !< Number of processes in current communicator
        integer :: west_rank                !< Previous rank ID in current communicator
        integer :: east_rank                !< Next rank ID in current communicator
        integer :: mpi_comm                 !< Current communicator
    end type cart_comm_1d

    type(cart_comm_1d), public :: comm_1d_x     !< Subcommunicator information in x-direction
    type(cart_comm_1d), public :: comm_1d_y     !< Subcommunicator information in y-direction
    type(cart_comm_1d), public :: comm_1d_z     !< Subcommunicator information in z-direction

    private

    public  :: mpi_topology_make
    public  :: mpi_topology_clean

    contains

    !>
    !> @brief       Destroy the communicator for cartesian topology.
    !>
    subroutine mpi_topology_clean()

        implicit none
        integer :: ierr

        call MPI_Comm_free(mpi_world_cart, ierr)

    end subroutine mpi_topology_clean

    !>
    !> @brief       Create the cartesian topology for the MPI processes and subcommunicators.
    !>
    subroutine mpi_topology_make()
        implicit none
        logical :: remain(0:2)
        integer :: ierr

        ! Create the cartesian topology.
        call MPI_Cart_create( MPI_COMM_WORLD    &!  input  | integer      | Input communicator (handle).
                            , 3                 &!  input  | integer      | Number of dimensions of Cartesian grid (integer).
                            , np_dim            &!  input  | integer(1:3) | Integer array of size ndims specifying the number of processes in each dimension.
                            , period            &!  input  | logical(1:3) | Logical array of size ndims specifying whether the grid is periodic (true=1) or not (false=0) in each dimension.
                            , .false.           &!  input  | logical      | Ranking may be reordered (true=1) or not (false=0) (logical).
                            , mpi_world_cart    &! *output | integer      | Communicator with new Cartesian topology (handle).
                            , ierr              &!  output | integer      | Fortran only: Error status
                            )

        ! Create subcommunicators and assign two neighboring processes in the x-direction.
        remain(0) = .true.
        remain(1) = .false.
        remain(2) = .false.
        call MPI_Cart_sub( mpi_world_cart, remain, comm_1d_x%mpi_comm, ierr)
        call MPI_Comm_rank(comm_1d_x%mpi_comm, comm_1d_x%myrank, ierr)
        call MPI_Comm_size(comm_1d_x%mpi_comm, comm_1d_x%nprocs, ierr)
        call MPI_Cart_shift(comm_1d_x%mpi_comm, 0, 1, comm_1d_x%west_rank, comm_1d_x%east_rank, ierr)

        ! Create subcommunicators and assign two neighboring processes in the y-direction
        remain(0) = .false.
        remain(1) = .true.
        remain(2) = .false.
        call MPI_Cart_sub( mpi_world_cart, remain, comm_1d_y%mpi_comm, ierr)
        call MPI_Comm_rank(comm_1d_y%mpi_comm, comm_1d_y%myrank, ierr)
        call MPI_Comm_size(comm_1d_y%mpi_comm, comm_1d_y%nprocs, ierr)
        call MPI_Cart_shift(comm_1d_y%mpi_comm, 0, 1, comm_1d_y%west_rank, comm_1d_y%east_rank, ierr)

        ! Create subcommunicators and assign two neighboring processes in the z-direction
        remain(0) = .false.
        remain(1) = .false.
        remain(2) = .true.
        call MPI_Cart_sub( mpi_world_cart, remain, comm_1d_z%mpi_comm, ierr)
        call MPI_Comm_rank(comm_1d_z%mpi_comm, comm_1d_z%myrank, ierr)
        call MPI_Comm_size(comm_1d_z%mpi_comm, comm_1d_z%nprocs, ierr)
        call MPI_Cart_shift(comm_1d_z%mpi_comm, 0, 1, comm_1d_z%west_rank, comm_1d_z%east_rank, ierr)

    end subroutine mpi_topology_make

end module mpi_topology
