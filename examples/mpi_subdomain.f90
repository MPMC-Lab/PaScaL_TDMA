!!======================================================================================================================
!> @file        mpi_subdomain.f90
!> @brief       This file contains a module of subdomains for the example problem of PaScaL_TDMA.
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
!> @brief       Module for building subdomains from the physical domain.
!> @details     This module has simulation parameters for subdomains and communication between the subdomains.
!>
module mpi_subdomain

    use MPI
    use global
    use mpi_topology, only : cart_comm_1d

    implicit none
    
    integer :: ierr
    !> @{ Grid numbers in the subdomain
    integer, public :: nx_sub,ny_sub,nz_sub
    !> @}
    !> @{ Grid indices of the assigned range
    integer, public :: ista, iend, jsta, jend, ksta, kend
    !> @}

    !> @{ Coordinates of grid points in the subdomain
    double precision, allocatable, dimension(:), public     :: x_sub, y_sub, Z_sub
    !> @}
    !> @{ Grid lengths in the subdomain
    double precision, allocatable, dimension(:), public     :: dmx_sub, dmy_sub, dmz_sub
    !> @}
    double precision, allocatable, dimension(:,:), public   :: thetaBC3_sub     !< B.C. of lower wall
    double precision, allocatable, dimension(:,:), public   :: thetaBC4_sub     !< B.C. of upper wall

    integer, allocatable, dimension(:), public  :: jmbc_index       !< Flag whether lower grid is empty(0) or not
    integer, allocatable, dimension(:), public  :: jpbc_index       !< Flag whether upper grid is empty(0) or not

    !> @{ Derived datatype for communication between x-neighbor subdomains
    integer :: ddtype_sendto_E, ddtype_recvfrom_W, ddtype_sendto_W, ddtype_recvfrom_E
    !> @}
    !> @{ Derived datatype for communication between y-neighbor subdomains
    integer :: ddtype_sendto_N, ddtype_recvfrom_S, ddtype_sendto_S, ddtype_recvfrom_N
    !> @}
    !> @{ Derived datatype for communication between z-neighbor subdomains
    integer :: ddtype_sendto_F, ddtype_recvfrom_B, ddtype_sendto_B, ddtype_recvfrom_F
    !> @}
    private

    public  :: mpi_subdomain_make
    public  :: mpi_subdomain_clean
    public  :: mpi_subdomain_make_ghostcell_ddtype
    public  :: mpi_subdomain_ghostcell_update
    public  :: mpi_subdomain_indices
    public  :: mpi_subdomain_mesh
    public  :: mpi_subdomain_initialization
    public  :: mpi_subdomain_boundary

    contains

    !>
    !> @brief       Prepare the subdomain and determine the size of the subdomain.
    !> @param       nprocs_in_x     Number of MPI processes in x-direction
    !> @param       myrank_in_x     Rank ID in x-direction
    !> @param       nprocs_in_y     Number of MPI processes in y-direction
    !> @param       myrank_in_y     Rank ID in y-direction
    !> @param       nprocs_in_z     Number of MPI processes in z-direction
    !> @param       myrank_in_z     Rank ID in z-direction
    !>
    subroutine mpi_subdomain_make(nprocs_in_x, myrank_in_x, &
                                         nprocs_in_y, myrank_in_y, &
                                         nprocs_in_z, myrank_in_z)
        implicit none
        integer, intent(in) :: nprocs_in_x, myrank_in_x, nprocs_in_y, myrank_in_y, nprocs_in_z, myrank_in_z

        ! Assigning grid numbers and grid indices of my subdomain.
        call para_range(1, nx-1, nprocs_in_x, myrank_in_x, ista, iend)
        nx_sub = iend - ista + 2
        call para_range(1, ny-1, nprocs_in_y, myrank_in_y, jsta, jend)
        ny_sub = jend - jsta + 2
        call para_range(1, nz-1, nprocs_in_z, myrank_in_z, ksta, kend)
        nz_sub = kend - ksta + 2
        
        ! Allocate subdomain variables.
        allocate( x_sub(0:nx_sub), dmx_sub(0:nx_sub))
        allocate( y_sub(0:ny_sub), dmy_sub(0:ny_sub))
        allocate( z_sub(0:nz_sub), dmz_sub(0:nz_sub))
        
        allocate( thetaBC3_sub(0:nx_sub, 0:nz_sub), thetaBC4_sub(0:nx_sub, 0:nz_sub))

        allocate(jmbc_index(0:ny_sub),jpbc_index(0:ny_sub))

    end subroutine mpi_subdomain_make

    !>
    !> @brief       Deallocate subdomain variables
    !>
    subroutine mpi_subdomain_clean
        implicit none

        deallocate(x_sub, dmx_sub)
        deallocate(y_sub, dmy_sub)
        deallocate(z_sub, dmz_sub)
        deallocate(thetaBC3_sub,thetaBC4_sub)

        deallocate(jmbc_index,jpbc_index)

    end subroutine mpi_subdomain_clean

    !>
    !> @brief       Build derived datatypes for subdomain communication using ghostcells.
    !>
    subroutine mpi_subdomain_make_ghostcell_ddtype

        implicit none
        integer :: sizes(0:2), subsizes(0:2), starts(0:2), ierr     ! Local variables for MPI_Type_create_subarray

        ! ddtype sending data to east MPI process (x+ neighbor)
        sizes    = (/nx_sub+1,ny_sub+1,nz_sub+1/)
        subsizes = (/      1,ny_sub+1,nz_sub+1/)
        starts   = (/nx_sub-1,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_sendto_E, ierr)
        call MPI_Type_commit(ddtype_sendto_E,ierr)

        ! ddtype receiving data from west MPI process (x- neighbor)
        starts   = (/      0,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_recvfrom_W, ierr)
        call MPI_Type_commit(ddtype_recvfrom_W,ierr)

        ! ddtype sending data to west MPI process (x- neighbor)
        starts   = (/      1,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_sendto_W, ierr)
        call MPI_Type_commit(ddtype_sendto_W,ierr)

        ! ddtype receiving data from east MPI process (x+ neighbor)
        starts   = (/  nx_sub,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_recvfrom_E, ierr)
        call MPI_Type_commit(ddtype_recvfrom_E,ierr)

        ! ddtype sending data to north MPI process (y+ neighbor)
        sizes    = (/nx_sub+1,ny_sub+1,nz_sub+1/)
        subsizes = (/nx_sub+1,      1,nz_sub+1/)
        starts   = (/      0,ny_sub-1,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_sendto_N, ierr)
        call MPI_Type_commit(ddtype_sendto_N,ierr)

        ! ddtype receiving data from south MPI process (y- neighbor)
        starts   = (/      0,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_recvfrom_S, ierr)
        call MPI_Type_commit(ddtype_recvfrom_S,ierr)

        ! ddtype sending data to south MPI process (y- neighbor)
        starts   = (/      0,      1,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_sendto_S, ierr)
        call MPI_Type_commit(ddtype_sendto_S,ierr)

        ! ddtype receiving data from north MPI process (y+ neighbor)
        starts   = (/      0,  ny_sub,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_recvfrom_N, ierr)
        call MPI_Type_commit(ddtype_recvfrom_N,ierr)

        ! ddtype sending data to forth MPI process (z+ neighbor)
        sizes    = (/nx_sub+1,ny_sub+1,nz_sub+1/)
        subsizes = (/nx_sub+1,ny_sub+1,      1/)
        starts   = (/      0,      0,nz_sub-1/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_sendto_F, ierr)
        call MPI_Type_commit(ddtype_sendto_F,ierr)

        ! ddtype receiving data from back MPI process (z- neighbor)
        starts   = (/      0,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_recvfrom_B, ierr)
        call MPI_Type_commit(ddtype_recvfrom_B,ierr)

        ! ddtype sending data to back MPI process (z- neighbor)
        starts   = (/      0,      0,      1/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_sendto_B, ierr)
        call MPI_Type_commit(  ddtype_sendto_B,ierr)

        ! ddtype receiving data from forth MPI process (z+ neighbor)
        starts   = (/      0,      0,  nz_sub/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_recvfrom_F, ierr)
        call MPI_Type_commit(ddtype_recvfrom_F,ierr)

    end subroutine mpi_subdomain_make_ghostcell_ddtype

    !>
    !> @brief       Update the values of boundary ghostcells through communication in all directions.
    !> @param       theta_sub       Variables to be updated
    !> @param       comm_1d_x       Subcommunicator in the x-direction
    !> @param       comm_1d_y       Subcommunicator in the y-direction
    !> @param       comm_1d_z       Subcommunicator in the z-direction
    !>
    subroutine mpi_subdomain_ghostcell_update(theta_sub, comm_1d_x, comm_1d_y, comm_1d_z)

        implicit none

        type(cart_comm_1d), intent(in)  :: comm_1d_x, comm_1d_y, comm_1d_z
        double precision, dimension(0:nx_sub, 0:ny_sub, 0:nz_sub), intent(inout) :: theta_sub

        integer :: ierr
        integer :: request(12)

        ! Update the ghostcells in the x-direction using derived datatypes and subcommunicator.
        call MPI_Isend(theta_sub,1, ddtype_sendto_E  , comm_1d_x%east_rank, 111, comm_1d_x%mpi_comm, request(1), ierr)
        call MPI_Irecv(theta_sub,1, ddtype_recvfrom_W, comm_1d_x%west_rank, 111, comm_1d_x%mpi_comm, request(2), ierr)
        call MPI_Isend(theta_sub,1, ddtype_sendto_W  , comm_1d_x%west_rank, 222, comm_1d_x%mpi_comm, request(3), ierr)
        call MPI_Irecv(theta_sub,1, ddtype_recvfrom_E, comm_1d_x%east_rank, 222, comm_1d_x%mpi_comm, request(4), ierr)
        
        ! Update the ghostcells in the y-direction using derived datatypes and subcommunicator.
        call MPI_Isend(theta_sub,1, ddtype_sendto_N  , comm_1d_y%east_rank, 111, comm_1d_y%mpi_comm, request(5), ierr)
        call MPI_Irecv(theta_sub,1, ddtype_recvfrom_S, comm_1d_y%west_rank, 111, comm_1d_y%mpi_comm, request(6), ierr)
        call MPI_Isend(theta_sub,1, ddtype_sendto_S  , comm_1d_y%west_rank, 222, comm_1d_y%mpi_comm, request(7), ierr)
        call MPI_Irecv(theta_sub,1, ddtype_recvfrom_N, comm_1d_y%east_rank, 222, comm_1d_y%mpi_comm, request(8), ierr)
        
        ! Update the ghostcells in the z-direction using derived datatypes and subcommunicator.
        call MPI_Isend(theta_sub,1, ddtype_sendto_F  , comm_1d_z%east_rank, 111, comm_1d_z%mpi_comm, request(9) , ierr)
        call MPI_Irecv(theta_sub,1, ddtype_recvfrom_B, comm_1d_z%west_rank, 111, comm_1d_z%mpi_comm, request(10), ierr)
        call MPI_Isend(theta_sub,1, ddtype_sendto_B  , comm_1d_z%west_rank, 222, comm_1d_z%mpi_comm, request(11), ierr)
        call MPI_Irecv(theta_sub,1, ddtype_recvfrom_F, comm_1d_z%east_rank, 222, comm_1d_z%mpi_comm, request(12), ierr)
        
        call MPI_Waitall(12, request, MPI_STATUSES_IGNORE, ierr)
        
    end subroutine mpi_subdomain_ghostcell_update

    !>
    !> @brief       Determine whether the next grids are empty(0) or not(1) only in the y-direction.
    !> @param       nprocs_in_y     Number of MPI processes in y-direction
    !> @param       myrank_in_y     Rank ID in y-direction
    !>
    subroutine mpi_subdomain_indices(myrank_in_y, nprocs_in_y)

        implicit none
        integer, intent(in) :: myrank_in_y, nprocs_in_y
    
        ! All values are initialized with 1, meaning all grids are not empty and effective.
        jmbc_index=1; jpbc_index=1
    
        ! Set the first jmbc_index to 0, having it empty in case of lower boundary domain.
        if(myrank_in_y==0) jmbc_index(1)=0
        ! Set the last jmbc_index to 0, having it empty in case of upper boundary domain.
        if(myrank_in_y==nprocs_in_y-1) jpbc_index(ny_sub-1)=0
    
        return
    end subroutine mpi_subdomain_indices

    !>
    !> @brief       Assign grid coordinates and lengths of subdomains.
    !> @param       myrank_in_x     Rank ID in x-direction
    !> @param       myrank_in_y     Rank ID in y-direction
    !> @param       myrank_in_z     Rank ID in z-direction
    !> @param       nprocs_in_x     Number of MPI processes in x-direction
    !> @param       nprocs_in_y     Number of MPI processes in y-direction
    !> @param       nprocs_in_z     Number of MPI processes in z-direction
    !>
    subroutine mpi_subdomain_mesh(myrank_in_x, myrank_in_y, myrank_in_z,  &
                                         nprocs_in_x, nprocs_in_y, nprocs_in_z)

        implicit none
        integer, intent(in) :: myrank_in_x,myrank_in_y,myrank_in_z,nprocs_in_x,nprocs_in_y,nprocs_in_z
    
        integer :: i,j,k

        ! X-direction: x_sub is for coordinates and dmx_sub is for grid lengths.
        dx = lx/dble(nx-1)
        do i = 0, nx_sub
            x_sub(i) = dble(i-1+ista-1)*dx
            dmx_sub(i)=dx
        end do

        ! Y-direction: y_sub is for coordinates and dmy_sub is for grid lengths.
        dy = ly/dble(ny)
        do j = 0, ny_sub
            y_sub(j) = dble(j+jsta-1)*dy
            dmy_sub(j)=dy
        end do

        ! Z-direction: z_sub is for coordinates and dmz_sub is for grid lengths.
        dz = lz/dble(nz-1)
        do k = 0, nz_sub
            z_sub(k) = dble(k-1+ksta-1)*dz
            dmz_sub(k)=dz
        end do

        ! For boundary grids in the x-direction. The same grid length is used due to periodic boundary conditions.
        if(myrank_in_x == 0) then
            dmx_sub(0)= dx
        endif
        if(myrank_in_x == nprocs_in_x-1) then
            dmx_sub(nx_sub)=dx
        endif

        ! For boundary grids in the y-direction. Half grid length is used for lower and upper boundary grids.
        if(myrank_in_y == 0) then
            dmy_sub(0)=dy/2.0d0
        endif
        if(myrank_in_y == nprocs_in_y-1) then
            dmy_sub(ny_sub)=dy/2.0d0
        endif 
    
        ! For boundary grids in the z-direction. The same grid length is used due to periodic boundary conditions.
        if(myrank_in_z == 0) then
            dmz_sub(0)=dz
        endif
        if(myrank_in_z == nprocs_in_z-1) then
            dmz_sub(nz_sub)=dz
        endif
    
        return
    
    end subroutine mpi_subdomain_mesh

    !>
    !> @brief       Initialize the values of the main variable in a subdomain.
    !> @param       theta_sub       Main variable to be solved
    !> @param       myrank_in_y     Rank ID in y-direction
    !> @param       nprocs_in_y     Number of MPI processes in y-direction
    !>
    subroutine mpi_subdomain_initialization(theta_sub, myrank_in_y, nprocs_in_y)
        implicit none
        integer, intent(in) :: myrank_in_y, nprocs_in_y
        double precision, intent(inout) :: theta_sub(0:nx_sub, 0:ny_sub, 0:nz_sub)
    
        integer :: i,j,k

        ! Initialize the main variable with a sine function and linearly changed values between the wall boundaries.
        do k = 0, nz_sub
            do j = 0, ny_sub
                do i = 0, nx_sub
                    theta_sub(i,j,k) = (theta_cold-theta_hot)/ly*Y_sub(j)+theta_hot   &
                                    & + sin(2*2.d0*PI/lx *X_sub(i))          &
                                    & * sin(2*2.d0*PI/lz *Z_sub(k))          &
                                    & * sin(2*2.d0*PI/ly *Y_sub(j))
                end do
            end do
        end do
    
        ! Initialize the values of the upper and lower boundary grids.
        do k = 0, nz_sub
            do i = 0, nx_sub
                if(myrank_in_y==0) theta_sub(i,0,k)  = theta_hot
                if(myrank_in_y==nprocs_in_y-1) theta_sub(i,ny_sub,k) = theta_cold
            end do
        end do
    
        return    
    end subroutine mpi_subdomain_initialization
    
    !>
    !> @brief       Assign the values of boundary grids in the subdomain.
    !> @param       theta_sub       Main variable to be solved
    !> @param       myrank_in_y     Rank ID in y-direction
    !> @param       nprocs_in_y     Number of MPI processes in y-direction
    !>
    subroutine mpi_subdomain_boundary(theta_sub, myrank_in_y, nprocs_in_y)
        implicit none
        integer, intent(in) :: myrank_in_y, nprocs_in_y
        double precision, intent(in) :: theta_sub(0:nx_sub, 0:ny_sub, 0:nz_sub)
    
        integer :: i,k

        ! Normal subdomain: Assign values of ghostcells to the boundary variables.
		do k = 0, nz_sub
			do i = 0, nx_sub
				thetaBC3_sub(i,k) = theta_sub(i,    0,k)
				thetaBC4_sub(i,k) = theta_sub(i,ny_sub,k)
			end do
		end do  

		! Lower boundary subdomain: Assign lower boundary condition to the boundary variables.
        if(myrank_in_y==0) then
            do k = 0, nz_sub
                do i = 0, nx_sub
                    thetaBC3_sub(i,k) = theta_hot
                end do
            end do  
		endif

		! Upper boundary subdomain: Assign upper boundary condition to the boundary variables.
        if(myrank_in_y==nprocs_in_y-1) then
            do k = 0, nz_sub
                do i = 0, nx_sub
                    thetaBC4_sub(i,k) = theta_cold
                end do
            end do
		endif

        return    
    end subroutine mpi_subdomain_boundary
    
end module mpi_subdomain
