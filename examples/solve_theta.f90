!======================================================================================================================
!> @file        solve_theta.f90
!> @brief       This file contains a solver subroutine for the example problem of PaScaL_TDMA.
!> @details     The target example problem is the three-dimensional time-dependent heat conduction problem 
!>              in a unit cube domain applied with the boundary conditions of vertically constant temperature 
!>              and horizontally periodic boundaries.
!>
!> @author      
!>              - Ki-Ha Kim (k-kiha@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
!>              - Mingyu Yang (yang926@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
!>              - Ji-Hoon Kang (jhkang@kisti.re.kr), Korea Institute of Science and Technology Information
!>              - Jung-Il Choi (jic@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
!>
!> @date        May 2023
!> @version     2.0
!> @par         Copyright
!>              Copyright (c) 2019-2023 Ki-Ha Kim, Mingyu Yang and Jung-Il choi, Yonsei University and 
!>              Ji-Hoon Kang, Korea Institute of Science and Technology Information, All rights reserved.
!> @par         License     
!>              This project is release under the terms of the MIT License (see LICENSE file).
!======================================================================================================================

module solve_theta

	implicit none

#ifdef _CUDA

    private

    ! Communication buffer
    !> @param       sbuf**          Buffer for MPI_Send in each direction and upper/lower GC
    !> @param       rbuf**          Buffer for MPI_Recv in each direction and upper/lower GC
    double precision, allocatable, dimension(:,:), device   :: sbuf_x0(:,:), sbuf_x1(:,:), sbuf_y0(:,:), &
                                                               sbuf_y1(:,:), sbuf_z0(:,:), sbuf_z1(:,:)
    double precision, allocatable, dimension(:,:), device   :: rbuf_x0(:,:), rbuf_x1(:,:), rbuf_y0(:,:), &
                                                               rbuf_y1(:,:), rbuf_z0(:,:), rbuf_z1(:,:)

    !> @param   dmx_sub_d       Grid lengths in the subdomain (device)
    !> @param   dmy_sub_d       Grid lengths in the subdomain (device)
    !> @param   dmz_sub_d       Grid lengths in the subdomain (device)
    !> @param   jpbc_index_d    Flag whether upper grid is empty(0) or not (device)
    !> @param   jmbc_index_d    Flag whether lower grid is empty(0) or not (device)
    !> @param   thetaBC3_sub_d  B.C. of lower wall (device)
    !> @param   thetaBC4_sub_d  B.C. of upper wall (device)

    double precision, allocatable, dimension(:), device    :: dmz_sub_d, dmy_sub_d, dmx_sub_d
    double precision, allocatable, dimension(:,:), device  :: thetaBC3_sub_d, thetaBC4_sub_d
    integer, allocatable, dimension(:), device             :: jpbc_index_d, jmbc_index_d
                                                       
    public  :: solve_theta_plan_many_cuda

#else
	public  :: solve_theta_plan_many
	public  :: solve_theta_plan_single
#endif

contains

#ifndef _CUDA
	!>
	!> @brief       An example solver for a single tridiagonal system of equations using PaScaL_TDMA.
	!> @details     This subroutine is for a single tridiagonal system of equations.
	!>              It solves the three-dimensional time-dependent heat conduction problem using PaScaL_TDMA solver.
	!>              PaScaL_TDMA plans are created for a single tridiagonal system of equations and
	!>              a single tridiagonal systems is solved line-by-line.
	!> @param       theta       Main 3-D variable to be solved
	!>
	subroutine solve_theta_plan_single(theta)
		use mpi
		use global
		use mpi_subdomain
		use mpi_topology
		use PaScaL_TDMA
		implicit none
		double precision, dimension(0:nx_sub, 0:ny_sub, 0:nz_sub), intent(inout) :: theta
		
		integer :: myrank, ierr
		integer :: time_step        ! Current time step
		double precision :: t_curr  ! Current simulation time

		! Loop and index variables
		integer :: i,j,k
		integer :: ip, jp, kp
		integer :: im, jm, km
		integer :: jem, jep

		! Temporary variables for coefficient computations
		double precision :: dedx1, dedx2, dedy3, dedy4, dedz5, dedz6    ! Derivative terms
		double precision :: viscous_e1, viscous_e2, viscous_e3, viscous ! Viscous terms
		double precision :: ebc_down, ebc_up, ebc                       ! Boundary terms
		double precision :: eAPI, eAMI, eACI                            ! Diffusion treatment terms in x-direction
		double precision :: eAPJ, eAMJ, eACJ                            ! Diffusion treatment terms in y-direction
		double precision :: eAPK, eAMK, eACK                            ! Diffusion treatment terms in z-direction
		double precision :: eRHS                                        ! From eAPI to eACK

		double precision, allocatable, dimension(:, :, :) :: rhs                    ! r.h.s. array
		double precision, allocatable, dimension(:) :: ap_1d, am_1d, ac_1d, ad_1d   ! Coefficient of ridiagonal matrix

		type(ptdma_plan_single)      :: px_single, py_single, pz_single ! Plan for a single tridiagonal system of equations

		call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr)

		t_curr = tStart
		dt = dtstart

		do time_step = 1, Tmax

			t_curr = t_curr + dt
			if(myrank==0) write(*,*) '[Main] Current time step = ', time_step
		
			! Calculating r.h.s
			allocate( rhs(0:nx_sub, 0:ny_sub, 0:nz_sub))
			do k = 1, nz_sub-1
				kp = k+1
				km = k-1

				do j = 1, ny_sub-1
					jp = j + 1
					jm = j - 1
					jep = jpbc_index(j)
					jem = jmbc_index(j)

					do i = 1, nx_sub-1
						ip = i+1
						im = i-1

						! Diffusion term
						dedx1 = (  theta(i ,j ,k ) - theta(im,j ,k )  )/dmx_sub(i )
						dedx2 = (  theta(ip,j ,k ) - theta(i ,j ,k )  )/dmx_sub(ip)  
						dedy3 = (  theta(i ,j ,k ) - theta(i ,jm,k )  )/dmy_sub(j )
						dedy4 = (  theta(i ,jp,k ) - theta(i ,j ,k )  )/dmy_sub(jp)
						dedz5 = (  theta(i ,j ,k ) - theta(i ,j ,km)  )/dmz_sub(k )
						dedz6 = (  theta(i ,j ,kp) - theta(i ,j ,k )  )/dmz_sub(kp)

						viscous_e1 = 1.d0/dx*(dedx2 - dedx1)
						viscous_e2 = 1.d0/dy*(dedy4 - dedy3)
						viscous_e3 = 1.d0/dz*(dedz6 - dedz5)
						viscous = 0.5d0*Ct*(viscous_e1 + viscous_e2 + viscous_e3) 
						
						! Boundary treatment for y-direction only
						ebc_down = 0.5d0*Ct/dy/dmy_sub(j)*thetaBC3_sub(i,k)
						ebc_up = 0.5d0*Ct/dy/dmy_sub(jp)*thetaBC4_sub(i,k)
						ebc = dble(1. - jem)*ebc_down + dble(1. - jep)*ebc_up

						! Diffusion term from incremental notation in next time step: x-direction
						eAPI = -0.5d0*Ct/dx/dmx_sub(ip)
						eAMI = -0.5d0*Ct/dx/dmx_sub(i )
						eACI =  0.5d0*Ct/dx*( 1.d0/dmx_sub(ip) + 1.d0/dmx_sub(i) )

						! Diffusion term from incremental notation in next time step: z-direction
						eAPK = -0.5d0*Ct/dz/dmz_sub(kp)
						eAMK = -0.5d0*Ct/dz/dmz_sub(k )
						eACK =  0.5d0*Ct/dz*( 1.d0/dmz_sub(kp) + 1.d0/dmz_sub(k) )

						! Diffusion term from incremental notation in next time step: y-direction
						eAPJ = -0.5d0*Ct/dy*( 1.d0/dmy_sub(jp) )*dble(jep)
						eAMJ = -0.5d0*Ct/dy*( 1.d0/dmy_sub(j ) )*dble(jem)
						eACJ =  0.5d0*Ct/dy*( 1.d0/dmy_sub(jp) + 1.d0/dmy_sub(j) )
						
						eRHS = eAPK*theta(i,j,kp) + eACK*theta(i,j,k) + eAMK*theta(i,j,km)      &
							& + eAPJ*theta(i,jp,k) + eACJ*theta(i,j,k) + eAMJ*theta(i,jm,k)      &
							& + eAPI*theta(ip,j,k) + eACI*theta(i,j,k) + eAMI*theta(im,j,k)

						! r.h.s.     
						rhs(i,j,k) = theta(i,j,k)/dt + viscous + ebc      &
								& - (theta(i,j,k)/dt + eRHS)
					enddo
				enddo
			enddo

			! Solve in the z-direction.
			allocate( ap_1d(1:nz_sub-1), am_1d(1:nz_sub-1), ac_1d(1:nz_sub-1), ad_1d(1:nz_sub-1) )

			! Create a PaScaL_TDMA plan for a single tridiagonal system.
			call PaScaL_TDMA_plan_single_create(pz_single, comm_1d_z%myrank, comm_1d_z%nprocs, comm_1d_z%mpi_comm, 0)

			! Build a coefficient matrix for a single tridiagonal system.
			do j = 1, ny_sub-1
				do i = 1, nx_sub-1
					do k = 1, nz_sub-1
						kp = k+1

						ap_1d(k) = -0.5d0*Ct/dz/dmz_sub(kp)*dt
						am_1d(k) = -0.5d0*Ct/dz/dmz_sub(k )*dt
						ac_1d(k) =  0.5d0*Ct/dz*( 1.d0/dmz_sub(kp) + 1.d0/dmz_sub(k) )*dt + 1.d0
						ad_1d(k) = rhs(i,j,k)*dt

					enddo
					! Solve a single tridiagonal system of equations under the defined plan with periodic boundary conditions.
					call PaScaL_TDMA_single_solve_cycle(pz_single, am_1d, ac_1d, ap_1d, ad_1d, nz_sub-1)
					! Return the solution to the r.h.s. line-by-line.
					rhs(i,j,1:nz_sub-1)=ad_1d(1:nz_sub-1)
				enddo
			enddo

			! Destroy the PaScaL_TDMA plan for a single tridiagonal system.
			call PaScaL_TDMA_plan_single_destroy(pz_single)
			deallocate( ap_1d, am_1d, ac_1d, ad_1d )

			! Solve in the y-direction.
			allocate( ap_1d(1:ny_sub-1), am_1d(1:ny_sub-1), ac_1d(1:ny_sub-1), ad_1d(1:ny_sub-1) )
			! Create a PaScaL_TDMA plan for a single tridiagonal system.
			call PaScaL_TDMA_plan_single_create(py_single, comm_1d_y%myrank, comm_1d_y%nprocs, comm_1d_y%mpi_comm, 0)

			! Build a coefficient matrix for a single tridiagonal system.
			do k = 1, nz_sub-1
				do i = 1, nx_sub-1
					do j = 1, ny_sub-1
						jp = j + 1
						jm = j - 1
						jep = jpbc_index(j)
						jem = jmbc_index(j)

						ap_1d(j) = -0.5d0*Ct/dy/dmy_sub(jp)*dble(jep)*dt
						am_1d(j) = -0.5d0*Ct/dy/dmy_sub(j )*dble(jem)*dt
						ac_1d(j) =  0.5d0*Ct/dy*( 1.d0/dmy_sub(jp) + 1.d0/dmy_sub(j) )*dt + 1.d0
						ad_1d(j) = rhs(i,j,k)
					end do
					! Solve a single tridiagonal system of equations under the defined plan.
					call PaScaL_TDMA_single_solve(py_single, am_1d, ac_1d, ap_1d, ad_1d, ny_sub-1)
					! Return the solution to r.h.s. line-by-line.
					rhs(i,1:ny_sub-1,k)=ad_1d(1:ny_sub-1)
				end do
			end do
			! Destroy the PaScaL_TDMA plan for a single tridiagonal system.
			call PaScaL_TDMA_plan_single_destroy(py_single)
			deallocate( ap_1d, am_1d, ac_1d, ad_1d )

			! Solve in the x-direction.
			allocate( ap_1d(1:nx_sub-1), am_1d(1:nx_sub-1), ac_1d(1:nx_sub-1), ad_1d(1:nx_sub-1) )
			! Create a PaScaL_TDMA plan for a single tridiagonal system.
			call PaScaL_TDMA_plan_single_create(px_single, comm_1d_x%myrank, comm_1d_x%nprocs, comm_1d_x%mpi_comm, 0)

			! Build a coefficient matrix for a single tridiagonal system.
			do k = 1, nz_sub-1
				do j = 1, ny_sub-1
					do i = 1, nx_sub-1
						ip = i+1
						im = i-1

						ap_1d(i) = -0.5d0*Ct/dx/dmx_sub(ip)*dt
						am_1d(i) = -0.5d0*Ct/dx/dmx_sub(i )*dt
						ac_1d(i) =  0.5d0*Ct/dx*( 1.d0/dmx_sub(ip) + 1.d0/dmx_sub(i) )*dt + 1.d0
						ad_1d(i) = rhs(i,j,k)
					enddo
					! Solve the single tridiagonal system of equations under the defined plan with periodic boundary conditions.
					call PaScaL_TDMA_single_solve_cycle(px_single, am_1d, ac_1d, ap_1d, ad_1d, nx_sub-1)
					! Return the solution to theta line-by-line.
					theta(1:nx_sub-1,j,k) = theta(1:nx_sub-1,j,k) + ad_1d(1:nx_sub-1)
				enddo
			enddo
			! Destroy the PaScaL_TDMA plan for a single tridiagonal system.
			call PaScaL_TDMA_plan_single_destroy(px_single)
			deallocate( ap_1d, am_1d, ac_1d, ad_1d )

			deallocate( rhs)

			! Update ghostcells from the solutions.
			call mpi_subdomain_ghostcell_update(theta, comm_1d_x, comm_1d_y, comm_1d_z)
		enddo

	end subroutine solve_theta_plan_single

	!>
	!> @brief       An example solver for many tridiagonal systems of equations using PaScaL_TDMA.
	!> @details     This subroutine is for many tridiagonal systems of equations.
	!>              It solves the the three-dimensional time-dependent heat conduction problem using PaScaL_TDMA.
	!>              PaScaL_TDMA plans are created for many tridiagonal systems of equations and
	!>              the many tridiagonal systems are solved plane-by-plane.
	!> @param       theta       Main 3-D variable to be solved
	!>
	subroutine solve_theta_plan_many(theta)

		use mpi
		use global
		use mpi_subdomain
		use mpi_topology
		use PaScaL_TDMA

		implicit none

		double precision, dimension(0:nx_sub, 0:ny_sub, 0:nz_sub), intent(inout) :: theta
		
		integer :: myrank, ierr
		integer :: time_step        ! Current time step
		double precision :: t_curr  ! Current simulation time

		! Loop and index variables
		integer :: i,j,k
		integer :: ip, jp, kp
		integer :: im, jm, km
		integer :: jem, jep

		! Temporary variables for coefficient computations
		double precision :: dedx1, dedx2, dedy3, dedy4, dedz5, dedz6    ! Derivative terms
		double precision :: viscous_e1, viscous_e2, viscous_e3, viscous ! Viscous terms
		double precision :: ebc_down, ebc_up, ebc                       ! Boundary terms
		double precision :: eAPI, eAMI, eACI                            ! Diffusion treatment terms in x-direction
		double precision :: eAPJ, eAMJ, eACJ                            ! Diffusion treatment terms in y-direction
		double precision :: eAPK, eAMK, eACK                            ! Diffusion treatment terms in z-direction
		double precision :: eRHS                                        ! From eAPI to eACK

		double precision, allocatable, dimension(:, :, :) :: rhs            ! r.h.s. array
		double precision, allocatable, dimension(:, :) :: ap, am, ac, ad    ! Coefficient of ridiagonal matrix

		type(ptdma_plan_many)     :: px_many, py_many, pz_many          ! Plan for many tridiagonal systems of equations

		call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr)

		t_curr = tStart
		dt = dtstart

		do time_step = 1, Tmax

			t_curr = t_curr + dt
			if(myrank==0) write(*,*) '[Main] Current time step = ', time_step
		
			! Calculating r.h.s.
			allocate( rhs(0:nx_sub, 0:ny_sub, 0:nz_sub))
			do k = 1, nz_sub-1
				kp = k+1
				km = k-1

				do j = 1, ny_sub-1
					jp = j + 1
					jm = j - 1
					jep = jpbc_index(j)
					jem = jmbc_index(j)

					do i = 1, nx_sub-1
						ip = i+1
						im = i-1

						! Diffusion term
						dedx1 = (  theta(i ,j ,k ) - theta(im,j ,k )  )/dmx_sub(i )
						dedx2 = (  theta(ip,j ,k ) - theta(i ,j ,k )  )/dmx_sub(ip)  
						dedy3 = (  theta(i ,j ,k ) - theta(i ,jm,k )  )/dmy_sub(j )
						dedy4 = (  theta(i ,jp,k ) - theta(i ,j ,k )  )/dmy_sub(jp)
						dedz5 = (  theta(i ,j ,k ) - theta(i ,j ,km)  )/dmz_sub(k )
						dedz6 = (  theta(i ,j ,kp) - theta(i ,j ,k )  )/dmz_sub(kp)

						viscous_e1 = 1.d0/dx*(dedx2 - dedx1)
						viscous_e2 = 1.d0/dy*(dedy4 - dedy3)
						viscous_e3 = 1.d0/dz*(dedz6 - dedz5)
						viscous = 0.5d0*Ct*(viscous_e1 + viscous_e2 + viscous_e3) 
						
						! Boundary treatment for the y-direction only
						ebc_down = 0.5d0*Ct/dy/dmy_sub(j)*thetaBC3_sub(i,k)
						ebc_up = 0.5d0*Ct/dy/dmy_sub(jp)*thetaBC4_sub(i,k)
						ebc = dble(1. - jem)*ebc_down + dble(1. - jep)*ebc_up

						! Diffusion term from incremental notation in next time step: x-direction
						eAPI = -0.5d0*Ct/dx/dmx_sub(ip)
						eAMI = -0.5d0*Ct/dx/dmx_sub(i )
						eACI =  0.5d0*Ct/dx*( 1.d0/dmx_sub(ip) + 1.d0/dmx_sub(i) )

						! Diffusion term from incremental notation in next time step: z-direction
						eAPK = -0.5d0*Ct/dz/dmz_sub(kp)
						eAMK = -0.5d0*Ct/dz/dmz_sub(k )
						eACK =  0.5d0*Ct/dz*( 1.d0/dmz_sub(kp) + 1.d0/dmz_sub(k) )

						! Diffusion term from incremental notation in next time step: y-direction
						eAPJ = -0.5d0*Ct/dy*( 1.d0/dmy_sub(jp) )*dble(jep)
						eAMJ = -0.5d0*Ct/dy*( 1.d0/dmy_sub(j ) )*dble(jem)
						eACJ =  0.5d0*Ct/dy*( 1.d0/dmy_sub(jp) + 1.d0/dmy_sub(j) )
						
						eRHS = eAPK*theta(i,j,kp) + eACK*theta(i,j,k) + eAMK*theta(i,j,km)      &
							& + eAPJ*theta(i,jp,k) + eACJ*theta(i,j,k) + eAMJ*theta(i,jm,k)      &
							& + eAPI*theta(ip,j,k) + eACI*theta(i,j,k) + eAMI*theta(im,j,k)

						! r.h.s. term 
						rhs(i,j,k) = theta(i,j,k)/dt + viscous + ebc      &
								& - (theta(i,j,k)/dt + eRHS)
					enddo
				enddo
			enddo

			! solve in the z-direction.
			allocate( ap(1:nx_sub-1, 1:nz_sub-1), am(1:nx_sub-1, 1:nz_sub-1), ac(1:nx_sub-1, 1:nz_sub-1), ad(1:nx_sub-1, 1:nz_sub-1) )

			! Create a PaScaL_TDMA plan for the tridiagonal systems.
			call PaScaL_TDMA_plan_many_create(pz_many, nx_sub-1, comm_1d_z%myrank, comm_1d_z%nprocs, comm_1d_z%mpi_comm)

			! Build a coefficient matrix for the tridiagonal systems into a 2D array.
			do j = 1, ny_sub-1
				do k = 1, nz_sub-1
					kp = k+1
					do i = 1, nx_sub-1
						ap(i,k) = -0.5d0*Ct/dz/dmz_sub(kp)*dt
						am(i,k) = -0.5d0*Ct/dz/dmz_sub(k )*dt
						ac(i,k) =  0.5d0*Ct/dz*( 1.d0/dmz_sub(kp) + 1.d0/dmz_sub(k) )*dt + 1.d0
						
						ad(i,k) = rhs(i,j,k)*dt
					enddo
				enddo

				! Solve the tridiagonal systems under the defined plan with periodic boundary conditions.
				call PaScaL_TDMA_many_solve_cycle(pz_many, am, ac, ap, ad,nx_sub-1,nz_sub-1)

				! Return the solution to the r.h.s. plane-by-plane
				do k = 1, nz_sub-1
					rhs(1:nx_sub-1,j,k)=ad(1:nx_sub-1,k)
				enddo
			enddo

			! Destroy the PaScaL_TDMA plan for the tridiagonal systems.
			call PaScaL_TDMA_plan_many_destroy(pz_many,comm_1d_z%nprocs)
			deallocate( ap, am, ac, ad )

			! solve in the x-direction.
			allocate( ap(1:nx_sub-1, 1:ny_sub-1), am(1:nx_sub-1, 1:ny_sub-1), ac(1:nx_sub-1, 1:ny_sub-1), ad(1:nx_sub-1, 1:ny_sub-1) )

			! Create a PaScaL_TDMA plan for the tridiagonal systems.
			call PaScaL_TDMA_plan_many_create(py_many, nx_sub-1, comm_1d_y%myrank, comm_1d_y%nprocs, comm_1d_y%mpi_comm)

			! Build a coefficient matrix for the tridiagonal systems into a 2D array.
			do k = 1, nz_sub-1
				do j = 1, ny_sub-1
					jp = j + 1
					jm = j - 1
					jep = jpbc_index(j)
					jem = jmbc_index(j)
					
					do i = 1, nx_sub-1
						ap(i,j) = -0.5d0*Ct/dy/dmy_sub(jp)*dble(jep)*dt
						am(i,j) = -0.5d0*Ct/dy/dmy_sub(j )*dble(jem)*dt
						ac(i,j) =  0.5d0*Ct/dy*( 1.d0/dmy_sub(jp) + 1.d0/dmy_sub(j) )*dt + 1.d0
						ad(i,j) = rhs(i,j,k)
					end do
				end do

				! Solve the tridiagonal systems under the defined plan.
				call PaScaL_TDMA_many_solve(py_many, am, ac, ap, ad, nx_sub-1, ny_sub-1)

				! Return the solution to the r.h.s. plane-by-plane.
				do j = 1, ny_sub-1
					rhs(1:nx_sub-1,j,k)=ad(1:nx_sub-1,j)
				enddo
			end do
			call PaScaL_TDMA_plan_many_destroy(py_many,comm_1d_y%nprocs)
			deallocate( ap, am, ac, ad )

			! solve in the y-direction.
			allocate( ap(1:ny_sub-1, 1:nx_sub-1), am(1:ny_sub-1, 1:nx_sub-1), ac(1:ny_sub-1, 1:nx_sub-1), ad(1:ny_sub-1, 1:nx_sub-1) )

			! Create a PaScaL_TDMA plan for the tridiagonal systems.
			call PaScaL_TDMA_plan_many_create(px_many, ny_sub-1, comm_1d_x%myrank, comm_1d_x%nprocs, comm_1d_x%mpi_comm)

			! Build a coefficient matrix for the tridiagonal systems into a 2D array.
			do k = 1, nz_sub-1
				do j = 1, ny_sub-1
					do i = 1, nx_sub-1
						ip = i+1
						im = i-1

						ap(j,i) = -0.5d0*Ct/dx/dmx_sub(ip)*dt
						am(j,i) = -0.5d0*Ct/dx/dmx_sub(i )*dt
						ac(j,i) =  0.5d0*Ct/dx*( 1.d0/dmx_sub(ip) + 1.d0/dmx_sub(i) )*dt + 1.d0
						ad(j,i) = rhs(i,j,k)
					enddo
				enddo
				! Solve the tridiagonal systems under the defined plan with periodic boundary conditions.
				call PaScaL_TDMA_many_solve_cycle(px_many, am, ac, ap, ad, ny_sub-1, nx_sub-1)

				! Return the solution to theta plane-by-plane.
				do j = 1, ny_sub-1
					theta(1:nx_sub-1,j,k) = theta(1:nx_sub-1,j,k) + ad(j,1:nx_sub-1)
				enddo

			enddo
			call PaScaL_TDMA_plan_many_destroy(px_many,comm_1d_x%nprocs)
			deallocate( ap, am, ac, ad )

			deallocate(rhs)

			! Update ghostcells from the solution.
			call mpi_subdomain_ghostcell_update(theta, comm_1d_x, comm_1d_y, comm_1d_z)
		end do

	end subroutine solve_theta_plan_many

#else
    !>
    !>
    !> @brief       An example solver for many tridiagonal systems of equations using PaScaL_TDMA with CUDA.
    !> @details     This subroutine is for many tridiagonal systems of equations.
    !>              It solves the the three-dimensional time-dependent heat conduction problem using PaScaL_TDMA.
    !>              PaScaL_TDMA plans are created for many tridiagonal systems of equations and
    !>              the many tridiagonal systems are solved plane-by-plane.
    !> @param       theta       Main 3-D variable to be solved
    !>
    subroutine solve_theta_plan_many_cuda(theta)

        use omp_lib
        use mpi
        use global, only : dx, dy, dz, dt, Ct, tStart, dtStart, Tmax
        use global, only : thread_in_x, thread_in_y, thread_in_z, thread_in_x_pascal, thread_in_y_pascal
        use mpi_subdomain, only : nx_sub, ny_sub, nz_sub
        use mpi_subdomain, only : dmx_sub, dmy_sub, dmz_sub
        use mpi_subdomain, only : jpbc_index, jmbc_index, thetaBC3_sub, thetaBC4_sub
        use mpi_topology, only : comm_1d_x, comm_1d_y, comm_1d_z
        use PaScaL_TDMA_cuda
        use cudafor
        use nvtx

        implicit none

        double precision, dimension(0:nx_sub, 0:ny_sub, 0:nz_sub), intent(inout) :: theta

        type(ptdma_plan_many_cuda)  :: px_many, pz_many, py_many  ! Plan for many tridiagonal systems of equations

        integer :: i, j, k
        integer :: myrank, ierr
        integer :: time_step        ! Current time step
        double precision :: t_curr  ! Current simulation time
    
        ! Temporary variables for coefficient computations
        double precision :: Ct_half_over_dx, Ct_half_over_dy, Ct_half_over_dz
        double precision, allocatable, dimension(:,:,:)         :: rhs              ! r.h.s. array
        double precision, allocatable, dimension(:,:,:), device :: theta_d, rhs_d   ! r.h.s. array

        double precision, allocatable, target, dimension(:), device :: ap_d, am_d, ac_d, ad_d   ! Coefficient of tridiagonal matrix
        double precision, pointer, dimension(:,:,:), device :: ap_ptr, am_ptr, ac_ptr, ad_ptr   ! Pointer to coefficient of tridiagonal matrix

        ! Block and thread dimension
        integer             :: block_in_x, block_in_y, block_in_z
        type(dim3)          :: blocks, threads, threads_in_pascal

        call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr)
    
        Ct_half_over_dx = 0.5d0*Ct/dx
        Ct_half_over_dy = 0.5d0*Ct/dy
        Ct_half_over_dz = 0.5d0*Ct/dz

        allocate( dmz_sub_d(0:nz_sub), dmy_sub_d(0:ny_sub), dmx_sub_d(0:nx_sub) )
        allocate( thetaBC3_sub_d(0:nx_sub, 0:nz_sub), thetaBC4_sub_d(0:nx_sub, 0:nz_sub) )
        allocate( jpbc_index_d(0:ny_sub), jmbc_index_d(0:ny_sub) )
    
        dmx_sub_d = dmx_sub
        dmy_sub_d = dmy_sub
        dmz_sub_d = dmz_sub
        jpbc_index_d = jpbc_index
        jmbc_index_d = jmbc_index
        thetaBC3_sub_d = thetaBC3_sub
        thetaBC4_sub_d = thetaBC4_sub

        ! Calculating r.h.s.
        allocate(   rhs_d(1:nx_sub-1, 1:ny_sub-1, 1:nz_sub-1));   rhs_d = 0.0d0
        allocate( theta_d(0:nx_sub  , 0:ny_sub  , 0:nz_sub  )); theta_d = 0.0d0

        allocate( sbuf_x0(0:ny_sub,0:nz_sub), sbuf_x1(0:ny_sub,0:nz_sub) ); sbuf_x0 = 0.0d0; sbuf_x1 = 0.0d0
        allocate( sbuf_y0(0:nx_sub,0:nz_sub), sbuf_y1(0:nx_sub,0:nz_sub) ); sbuf_y0 = 0.0d0; sbuf_y1 = 0.0d0
        allocate( sbuf_z0(0:nx_sub,0:ny_sub), sbuf_z1(0:nx_sub,0:ny_sub) ); sbuf_z0 = 0.0d0; sbuf_z1 = 0.0d0

        allocate( rbuf_x0(0:ny_sub,0:nz_sub), rbuf_x1(0:ny_sub,0:nz_sub) ); rbuf_x0 = 0.0d0; rbuf_x1 = 0.0d0
        allocate( rbuf_y0(0:nx_sub,0:nz_sub), rbuf_y1(0:nx_sub,0:nz_sub) ); rbuf_y0 = 0.0d0; rbuf_y1 = 0.0d0
        allocate( rbuf_z0(0:nx_sub,0:ny_sub), rbuf_z1(0:nx_sub,0:ny_sub) ); rbuf_z0 = 0.0d0; rbuf_z1 = 0.0d0

        ! Setting the thread and block
        block_in_x = (nx_sub-1)/thread_in_x
        if(block_in_x.eq.0 .or. mod((nx_sub-1), thread_in_x)) then
            print '(a,i5,a,i5)', '[Error] ny_sub-1 should be a multiple of thread_in_x. &
                                   thread_in_x = ',thread_in_x,', nx_sub-1 = ',nx_sub-1
            call MPI_Finalize(ierr)
            stop
        endif

        block_in_y = (ny_sub-1)/thread_in_y
        if(block_in_y.eq.0 .or. mod((ny_sub-1), thread_in_y)) then
            print '(a,i5,a,i5)', '[Error] nz_sub-1 should be a multiple of thread_in_y. &
                                   thread_in_y = ',thread_in_y,', ny_sub-1 = ',ny_sub-1
            call MPI_Finalize(ierr)
            stop
        endif

        block_in_z = (nz_sub-1)/thread_in_z
        if(block_in_z.eq.0 .or. mod((nz_sub-1), thread_in_z)) then
            print '(a,i5,a,i5)', '[Error] nz_sub-1 should be a multiple of thread_in_z. &
                                   thread_in_z = ',thread_in_z,', nz_sub-1 = ',nz_sub-1
            call MPI_Finalize(ierr)
            stop
        endif

        blocks  = dim3(block_in_x, block_in_y, block_in_z)
        threads = dim3(thread_in_x, thread_in_y, thread_in_z)

        threads_in_pascal = dim3(thread_in_x_pascal, thread_in_y_pascal, 1)
    
        ! Simulation begins
        t_curr = tStart
        dt = dtstart

        ! Create a PaScaL_TDMA plan for the tridiagonal systems.
        call PaScaL_TDMA_plan_many_create_cuda( pz_many, nx_sub-1, ny_sub-1, nz_sub-1, comm_1d_z%myrank, &
                                                comm_1d_z%nprocs, comm_1d_z%mpi_comm, threads_in_pascal)
        call PaScaL_TDMA_plan_many_create_cuda( py_many, nx_sub-1, nz_sub-1, ny_sub-1, comm_1d_y%myrank, &
                                                comm_1d_y%nprocs, comm_1d_y%mpi_comm, threads_in_pascal)
        call PaScaL_TDMA_plan_many_create_cuda( px_many, ny_sub-1, nz_sub-1, nx_sub-1, comm_1d_x%myrank, &
                                                comm_1d_x%nprocs, comm_1d_x%mpi_comm, threads_in_pascal)

        allocate( ap_d((nx_sub-1)*(ny_sub-1)*(nz_sub-1)) ); ap_d = 0.0d0
        allocate( ac_d((nx_sub-1)*(ny_sub-1)*(nz_sub-1)) ); ac_d = 0.0d0
        allocate( am_d((nx_sub-1)*(ny_sub-1)*(nz_sub-1)) ); am_d = 0.0d0
        allocate( ad_d((nx_sub-1)*(ny_sub-1)*(nz_sub-1)) ); ad_d = 0.0d0

        theta_d = theta

        do time_step = 1, Tmax

            t_curr = t_curr + dt
            if(myrank==0) write(*,*) '[Main] Current time step = ', time_step

            call build_RHS_cuda<<<blocks, threads>>>(theta_d, rhs_d, &
                                                    Ct_half_over_dx, Ct_half_over_dy, Ct_half_over_dz)

            ! solve in the z-direction.
            ap_ptr(1:nx_sub-1, 1:ny_sub-1, 1:nz_sub-1) => ap_d
            ac_ptr(1:nx_sub-1, 1:ny_sub-1, 1:nz_sub-1) => ac_d
            am_ptr(1:nx_sub-1, 1:ny_sub-1, 1:nz_sub-1) => am_d
            ad_ptr(1:nx_sub-1, 1:ny_sub-1, 1:nz_sub-1) => ad_d

            ! Build a coefficient matrix for the tridiagonal systems into a 2D array.
            call build_LHSz_cuda<<<blocks, threads>>>(ap_ptr, am_ptr, ac_ptr, ad_ptr, rhs_d, Ct_half_over_dz, dt)

            !Solve the tridiagonal systems under the defined plan with periodic boundary conditions.
            call PaScaL_TDMA_many_solve_cycle_cuda( pz_many, am_ptr, ac_ptr, ap_ptr, ad_ptr)

            ! Return the solution to the r.h.s. plane-by-plane
            call copy_ijk2ijk<<<blocks, threads>>>(ad_ptr, rhs_d)

            nullify( ap_ptr, ac_ptr, am_ptr, ad_ptr )
            ! solve in the y-direction.
            ap_ptr(1:nx_sub-1, 1:nz_sub-1, 1:ny_sub-1) => ap_d
            ac_ptr(1:nx_sub-1, 1:nz_sub-1, 1:ny_sub-1) => ac_d
            am_ptr(1:nx_sub-1, 1:nz_sub-1, 1:ny_sub-1) => am_d
            ad_ptr(1:nx_sub-1, 1:nz_sub-1, 1:ny_sub-1) => ad_d

            ! Build a coefficient matrix for the tridiagonal systems into a 2D array.
            call build_LHSy_cuda<<<blocks, threads>>>(ap_ptr, am_ptr, ac_ptr, ad_ptr, rhs_d, Ct_half_over_dy, dt)

            ! Solve the tridiagonal systems under the defined plan.
            call PaScaL_TDMA_many_solve_cuda(py_many, am_ptr, ac_ptr, ap_ptr, ad_ptr)
            
            ! Return the solution to the r.h.s. plane-by-plane.
            call transpose_ikj2ijk<<<blocks, threads>>>(ad_ptr, rhs_d)

            nullify( ap_ptr, ac_ptr, am_ptr, ad_ptr )

            ! solve in the x-direction.
            ap_ptr(1:ny_sub-1, 1:nz_sub-1, 1:nx_sub-1) => ap_d
            ac_ptr(1:ny_sub-1, 1:nz_sub-1, 1:nx_sub-1) => ac_d
            am_ptr(1:ny_sub-1, 1:nz_sub-1, 1:nx_sub-1) => am_d
            ad_ptr(1:ny_sub-1, 1:nz_sub-1, 1:nx_sub-1) => ad_d

            ! Build a coefficient matrix for the tridiagonal systems into a 2D array.
            call build_LHSx_cuda<<<blocks, threads>>>(ap_ptr, am_ptr, ac_ptr, ad_ptr, rhs_d, Ct_half_over_dx, dt)

            ! Solve the tridiagonal systems under the defined plan with periodic boundary conditions.
            call PaScaL_TDMA_many_solve_cycle_cuda( px_many, am_ptr, ac_ptr, ap_ptr, ad_ptr)

            ! Return the solution to theta plane-by-plane.
            call update_theta_cuda<<<blocks, threads>>>(ad_ptr, theta_d)

            nullify( ap_ptr, ac_ptr, am_ptr, ad_ptr )

            ! ! Update ghostcells from the solution.
            call ghostcell_update_cuda(theta_d)
        end do

        theta = theta_d

        ! Destroy the PaScaL_TDMA plan for the tridiagonal systems.
        call PaScaL_TDMA_plan_many_destroy_cuda(pz_many)
        call PaScaL_TDMA_plan_many_destroy_cuda(px_many)
        call PaScaL_TDMA_plan_many_destroy_cuda(py_many)

        deallocate( ap_d, am_d, ac_d, ad_d )
        deallocate( rhs_d, theta_d )

        deallocate( sbuf_x0, sbuf_x1 )
        deallocate( sbuf_y0, sbuf_y1 )
        deallocate( sbuf_z0, sbuf_z1 )
        deallocate( rbuf_x0, rbuf_x1 )
        deallocate( rbuf_y0, rbuf_y1 )
        deallocate( rbuf_z0, rbuf_z1 )

        deallocate( dmz_sub_d, dmy_sub_d, dmx_sub_d )
        deallocate( thetaBC3_sub_d, thetaBC4_sub_d )
        deallocate( jpbc_index_d, jmbc_index_d )
    

    end subroutine solve_theta_plan_many_cuda

    !>
    !>
    !> @brief       Ghost cell update in device
    !> @details     This subroutine is to update the ghost-cell in the device.
    !>              Allocated buffers in device memory are used. 
    !> @param       Value_sub_d     Main 3-D variable to be updated with GC
    !>
    subroutine ghostcell_update_cuda(Value_sub_d)

        use mpi
        use global
        use mpi_subdomain
        use mpi_topology
        use cudafor

        implicit none

        double precision, device, dimension(0:nx_sub, 0:ny_sub, 0:nz_sub), intent(inout)  :: Value_sub_d
        integer :: i, j, k
        integer :: ierr, request(4)
        integer :: nprocs
        integer :: istat(MPI_STATUS_SIZE)
                
        !X-Direction        
        !$cuf kernel do(2) <<< *,* >>>
        do k = 0, nz_sub
        do j = 0, ny_sub
            if(comm_1d_x%west_rank.ne.MPI_PROC_NULL) then
                sbuf_x0(j,k) = Value_sub_d(1       ,j,k)
            endif
            if(comm_1d_x%east_rank.ne.MPI_PROC_NULL) then
                sbuf_x1(j,k) = Value_sub_d(nx_sub-1,j,k)
            endif
        enddo
        enddo

        if( comm_1d_x%nprocs.eq.1 .and. period(0).eqv..true. ) then
            !$cuf kernel do(2) <<< *,* >>>
            do k = 0, nz_sub
            do j = 0, ny_sub
                rbuf_x1(j,k) = sbuf_x0(j,k)
                rbuf_x0(j,k) = sbuf_x1(j,k) 
            enddo
            enddo
        else
            ierr = cudaStreamSynchronize()
            call MPI_Isend( sbuf_x0, (ny_sub+1)*(nz_sub+1), MPI_DOUBLE_PRECISION, &
                            comm_1d_x%west_rank, 111, comm_1d_x%mpi_comm, request(1), ierr)
            call MPI_Irecv( rbuf_x1, (ny_sub+1)*(nz_sub+1), MPI_DOUBLE_PRECISION, &
                            comm_1d_x%east_rank, 111, comm_1d_x%mpi_comm, request(2), ierr)
            call MPI_Irecv( rbuf_x0, (ny_sub+1)*(nz_sub+1), MPI_DOUBLE_PRECISION, &
                            comm_1d_x%west_rank, 222, comm_1d_x%mpi_comm, request(3), ierr)
            call MPI_Isend( sbuf_x1, (ny_sub+1)*(nz_sub+1), MPI_DOUBLE_PRECISION, &
                            comm_1d_x%east_rank, 222, comm_1d_x%mpi_comm, request(4), ierr)
            call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)
        endif

        !$cuf kernel do(2) <<< *,* >>>
        do k = 0, nz_sub
        do j = 0, ny_sub
            if(comm_1d_x%west_rank.ne.MPI_PROC_NULL) then
                Value_sub_d(0     ,j,k) = rbuf_x0(j,k)
            endif
            if(comm_1d_x%east_rank.ne.MPI_PROC_NULL) then
                Value_sub_d(nx_sub,j,k) = rbuf_x1(j,k)
            endif
        enddo
        enddo

        !Y-Direction    
        if(comm_1d_y%west_rank.ne.MPI_PROC_NULL) then
            !$cuf kernel do(2) <<< *,* >>>
            do k = 0, nz_sub
            do i = 0, nx_sub
                sbuf_y0(i,k) = Value_sub_d(i,1       ,k)
            enddo
            enddo
        endif
        if(comm_1d_y%east_rank.ne.MPI_PROC_NULL) then
            !$cuf kernel do(2) <<< *,* >>>
            do k = 0, nz_sub
            do i = 0, nx_sub
                sbuf_y1(i,k) = Value_sub_d(i,ny_sub-1,k)
            enddo
            enddo
        endif

        if( comm_1d_y%nprocs.eq.1 .and. period(1).eqv..true. ) then 
            !$cuf kernel do(2) <<< *,* >>>
            do k = 0, nz_sub
            do i = 0, nx_sub
                rbuf_y1(i,k) = sbuf_y0(i,k)
                rbuf_y0(i,k) = sbuf_y1(i,k) 
            enddo
            enddo
        else
            ierr = cudaStreamSynchronize()
            call MPI_Isend( sbuf_y0, (nx_sub+1)*(nz_sub+1), MPI_DOUBLE_PRECISION, &
                            comm_1d_y%west_rank, 333, comm_1d_y%mpi_comm, request(1), ierr)
            call MPI_Irecv(rbuf_y1, (nx_sub+1)*(nz_sub+1), MPI_DOUBLE_PRECISION, &
                            comm_1d_y%east_rank, 333, comm_1d_y%mpi_comm, request(2), ierr)
            call MPI_Irecv(rbuf_y0, (nx_sub+1)*(nz_sub+1), MPI_DOUBLE_PRECISION, &
                            comm_1d_y%west_rank, 444, comm_1d_y%mpi_comm, request(3), ierr)
            call MPI_Isend( sbuf_y1, (nx_sub+1)*(nz_sub+1), MPI_DOUBLE_PRECISION, &
                            comm_1d_y%east_rank, 444, comm_1d_y%mpi_comm, request(4), ierr)
            call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)
        endif

        if(comm_1d_y%west_rank.ne.MPI_PROC_NULL) then
        !$cuf kernel do(2) <<< *,* >>>
            do k = 0, nz_sub
            do i = 0, nx_sub
                Value_sub_d(i,0     ,k) = rbuf_y0(i,k)
            enddo
            enddo
        endif
        if(comm_1d_y%east_rank.ne.MPI_PROC_NULL) then
        !$cuf kernel do(2) <<< *,* >>>
            do k = 0, nz_sub
            do i = 0, nx_sub
                Value_sub_d(i,ny_sub,k) = rbuf_y1(i,k)
            enddo
            enddo
        endif

        !Z-Direction
        if(comm_1d_z%west_rank.ne.MPI_PROC_NULL) then
            !$cuf kernel do(2) <<< *,* >>>
            do j = 0, ny_sub
            do i = 0, nx_sub
                sbuf_z0(i,j) = Value_sub_d(i,j,1       )
            enddo
            enddo
        endif
        if(comm_1d_z%east_rank.ne.MPI_PROC_NULL) then
            !$cuf kernel do(2) <<< *,* >>>
            do j = 0, ny_sub
            do i = 0, nx_sub
                sbuf_z1(i,j) = Value_sub_d(i,j,nz_sub-1)
            enddo
            enddo
        endif
       
        if( comm_1d_z%nprocs.eq.1 .and. period(2).eqv..true. ) then
            !$cuf kernel do(2) <<< *,* >>>
            do j = 0, ny_sub
            do i = 0, nx_sub
                rbuf_z1(i,j) = sbuf_z0(i,j)
                rbuf_z0(i,j) = sbuf_z1(i,j) 
            enddo
            enddo
        else
            ierr = cudaStreamSynchronize()
            call MPI_Isend( sbuf_z0, (nx_sub+1)*(ny_sub+1), MPI_DOUBLE_PRECISION, &
                            comm_1d_z%west_rank, 555, comm_1d_z%mpi_comm, request(1), ierr)
            call MPI_Irecv( rbuf_z1, (nx_sub+1)*(ny_sub+1), MPI_DOUBLE_PRECISION, &
                            comm_1d_z%east_rank, 555, comm_1d_z%mpi_comm, request(2), ierr)
            call MPI_Irecv( rbuf_z0, (nx_sub+1)*(ny_sub+1), MPI_DOUBLE_PRECISION, &
                            comm_1d_z%west_rank, 666, comm_1d_z%mpi_comm, request(3), ierr)
            call MPI_Isend( sbuf_z1, (nx_sub+1)*(ny_sub+1), MPI_DOUBLE_PRECISION, &
                            comm_1d_z%east_rank, 666, comm_1d_z%mpi_comm, request(4), ierr)
            call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)
        endif

        if(comm_1d_z%west_rank.ne.MPI_PROC_NULL) then
            !$cuf kernel do(2) <<< *,* >>>
            do j = 0, ny_sub
            do i = 0, nx_sub
                Value_sub_d(i,j,0     ) = rbuf_z0(i,j)
            enddo
            enddo
        endif
        if(comm_1d_z%east_rank.ne.MPI_PROC_NULL) then
            !$cuf kernel do(2) <<< *,* >>>
            do j = 0, ny_sub
            do i = 0, nx_sub
                Value_sub_d(i,j,nz_sub) = rbuf_z1(i,j)
            enddo
            enddo
        endif
        
    end subroutine ghostcell_update_cuda

    !>
    !> @brief   Build RHS term using initial conditions and geometry information
    !> @param   theda_d         Updated solution at the previous time step (device)
    !> @param   rhd_d           RHS term (device)
    !> @param   coefx           Ct_half_over_dx (device)
    !> @param   coefy           Ct_half_over_dy (device)
    !> @param   coefz           Ct_half_over_dz (device)
    !>
    attributes(global) subroutine build_RHS_cuda(theta_d, rhs_d, coefx, coefy, coefz)

        implicit none

        double precision, device, intent(in)    :: theta_d(0:, 0:, 0:)       ! r.h.s. array
        double precision, device, intent(inout) :: rhs_d(:, :, :)            ! r.h.s. array
        double precision, value, intent(in)     :: coefx, coefy, coefz

        ! Loop and index variables
        integer :: i,j,k
        integer :: jem, jep

        ! Temporary variables for coefficient computations
        double precision :: dedx1, dedx2, dedy3, dedy4, dedz5, dedz6    ! Derivative terms
        double precision :: viscous                                     ! Viscous terms
        double precision :: ebc_down, ebc_up, ebc                       ! Boundary terms
        double precision :: eAPI, eAMI, eACI                            ! Diffusion treatment terms in x-direction
        double precision :: eAPJ, eAMJ, eACJ                            ! Diffusion treatment terms in y-direction
        double precision :: eAPK, eAMK, eACK                            ! Diffusion treatment terms in z-direction
        double precision :: eRHS                                        ! From eAPI to eACK
        double precision :: inv_dmx_sub_i, inv_dmx_sub_ip, inv_dmy_sub_j, inv_dmy_sub_jp, inv_dmz_sub_k, inv_dmz_sub_kp
        double precision :: thetaijk, thetaip, thetaim, thetajp, thetajm, thetakp, thetakm
        double precision :: thetaBC3, thetaBC4

        i = (blockidx%x-1) * blockdim%x + threadidx%x
        j = (blockidx%y-1) * blockdim%y + threadidx%y
        k = (blockidx%z-1) * blockdim%z + threadidx%z

        thetaijk = theta_d(i,  j,  k  )
        thetaip  = theta_d(i+1,j,  k  )
        thetaim  = theta_d(i-1,j,  k  )
        thetajp  = theta_d(i,  j+1,k  )
        thetajm  = theta_d(i,  j-1,k  )
        thetakp  = theta_d(i,  j,  k+1)
        thetakm  = theta_d(i,  j,  k-1)

        thetaBC3 = thetaBC3_sub_d(i,k)
        thetaBC4 = thetaBC4_sub_d(i,k)

        inv_dmx_sub_i  = 1.0d0 / dmx_sub_d(i  )
        inv_dmx_sub_ip = 1.0d0 / dmx_sub_d(i+1)
        inv_dmy_sub_j  = 1.0d0 / dmy_sub_d(j  )
        inv_dmy_sub_jp = 1.0d0 / dmy_sub_d(j+1)
        inv_dmz_sub_k  = 1.0d0 / dmz_sub_d(k  )
        inv_dmz_sub_kp = 1.0d0 / dmz_sub_d(k+1)

        jep = jpbc_index_d(j)
        jem = jmbc_index_d(j)

        ! Diffusion term
        dedx1 = (  thetaijk - thetaim  )*inv_dmx_sub_i
        dedx2 = (  thetaip  - thetaijk )*inv_dmx_sub_ip
        dedy3 = (  thetaijk - thetajm  )*inv_dmy_sub_j
        dedy4 = (  thetajp  - thetaijk )*inv_dmy_sub_jp
        dedz5 = (  thetaijk - thetakm  )*inv_dmz_sub_k
        dedz6 = (  thetakp  - thetaijk )*inv_dmz_sub_kp

        viscous = coefx*(dedx2 - dedx1) + coefy*(dedy4 - dedy3) + coefz*(dedz6 - dedz5)
        
        ! Boundary treatment for the y-direction only
        ebc_down = coefy*inv_dmy_sub_j *thetaBC3
        ebc_up   = coefy*inv_dmy_sub_jp*thetaBC4
        ebc = dble(1-jem)*ebc_down + dble(1-jep)*ebc_up

        ! Diffusion term from incremental notation in next time step: x-direction
        eAPI = -coefx*inv_dmx_sub_ip
        eAMI = -coefx*inv_dmx_sub_i
        eACI =  coefx*( inv_dmx_sub_ip + inv_dmx_sub_i )

        ! Diffusion term from incremental notation in next time step: z-direction
        eAPK = -coefz*inv_dmz_sub_kp
        eAMK = -coefz*inv_dmz_sub_k
        eACK =  coefz*( inv_dmz_sub_kp + inv_dmz_sub_k )

        ! Diffusion term from incremental notation in next time step: y-direction
        eAPJ = -coefy*inv_dmy_sub_jp*dble(jep)
        eAMJ = -coefy*inv_dmy_sub_j *dble(jem)
        eACJ =  coefy*( inv_dmy_sub_jp + inv_dmy_sub_j )

        eRHS  = eAPK*thetakp + eACK*thetaijk + eAMK*thetakm      &
            & + eAPJ*thetajp + eACJ*thetaijk + eAMJ*thetajm      &
            & + eAPI*thetaip + eACI*thetaijk + eAMI*thetaim

        ! r.h.s. term 
        rhs_d(i,j,k) = viscous + ebc - eRHS

    end subroutine build_RHS_cuda

    !>
    !> @brief   Build tri-diagonal matrix for solving in the z-direction
    !> @param   ap_ptr          Upper diagonal term (device)
    !> @param   am_ptr          Lower diagonal term (device)
    !> @param   ac_ptr          Diagonal term (device)
    !> @param   ad_ptr          Updated RHS term (device)
    !> @param   rhs_d           RHS term (device)
    !> @param   coefz           Ct_half_over_dz (device)
    !> @param   dt              Time step (device)
    !>
    attributes(global) subroutine build_LHSz_cuda(ap_ptr, am_ptr, ac_ptr, ad_ptr, rhs_d, coefz, dt)

        implicit none

        double precision, device, intent(inout) :: ap_ptr(:, :, :), am_ptr(:, :, :)
        double precision, device, intent(inout) :: ac_ptr(:, :, :), ad_ptr(:, :, :)
        double precision, device, intent(in)    :: rhs_d(:, :, :)            ! r.h.s. array
        double precision, value, intent(in)     :: coefz, dt

        integer :: i,j,k
        double precision :: inv_dmz_sub_k, inv_dmz_sub_kp

        i = (blockidx%x-1) * blockdim%x + threadidx%x
        j = (blockidx%y-1) * blockdim%y + threadidx%y
        k = (blockidx%z-1) * blockdim%z + threadidx%z

        inv_dmz_sub_k  = 1.0d0 / dmz_sub_d(k  )
        inv_dmz_sub_kp = 1.0d0 / dmz_sub_d(k+1)

        ap_ptr(i,j,k) = -coefz*inv_dmz_sub_kp*dt
        am_ptr(i,j,k) = -coefz*inv_dmz_sub_k *dt
        ac_ptr(i,j,k) =  coefz*( inv_dmz_sub_kp + inv_dmz_sub_k )*dt + 1.d0
        ad_ptr(i,j,k) =  rhs_d(i,j,k) * dt

    end subroutine build_LHSz_cuda

    !>
    !> @brief   Build tri-diagonal matrix for solving in the y-direction
    !> @param   ap_ptr          Upper diagonal term (device)
    !> @param   am_ptr          Lower diagonal term (device)
    !> @param   ac_ptr          Diagonal term (device)
    !> @param   ad_ptr          Updated RHS term (device)
    !> @param   rhs_d           RHS term (device)
    !> @param   coefy           Ct_half_over_dy (device)
    !> @param   dt              Time step (device)
    !>
    attributes(global) subroutine build_LHSy_cuda(ap_ptr, am_ptr, ac_ptr, ad_ptr, rhs_d, coefy, dt)

        implicit none

        double precision, device, intent(inout) :: ap_ptr(:, :, :), am_ptr(:, :, :)
        double precision, device, intent(inout) :: ac_ptr(:, :, :), ad_ptr(:, :, :)
        double precision, device, intent(in)    :: rhs_d(:, :, :)            ! r.h.s. array
        double precision, value, intent(in)     :: coefy, dt

        integer :: i,j,k
        integer :: jep, jem
        double precision :: inv_dmy_sub_j, inv_dmy_sub_jp

        i = (blockidx%x-1) * blockdim%x + threadidx%x
        j = (blockidx%y-1) * blockdim%y + threadidx%y
        k = (blockidx%z-1) * blockdim%z + threadidx%z

        inv_dmy_sub_j  = 1.0d0 / dmy_sub_d(j  )
        inv_dmy_sub_jp = 1.0d0 / dmy_sub_d(j+1)

        jep = jpbc_index_d(j)
        jem = jmbc_index_d(j)
        
        ap_ptr(i,k,j) = -coefy*inv_dmy_sub_jp*dble(jep)*dt
        am_ptr(i,k,j) = -coefy*inv_dmy_sub_j *dble(jem)*dt
        ac_ptr(i,k,j) =  coefy*( inv_dmy_sub_jp + inv_dmy_sub_j )*dt + 1.0d0
        ad_ptr(i,k,j) =  rhs_d(i,j,k)

    end subroutine build_LHSy_cuda

    !>
    !> @brief   Build tri-diagonal matrix for solving in the z-direction
    !> @param   ap_ptr          Upper diagonal term (device)
    !> @param   am_ptr          Lower diagonal term (device)
    !> @param   ac_ptr          Diagonal term (device)
    !> @param   ad_ptr          Updated RHS term (device)
    !> @param   rhs_d           RHS term (device)
    !> @param   coefx           Ct_half_over_dx (device)
    !> @param   dt              Time step (device)
    !>
    attributes(global) subroutine build_LHSx_cuda(ap_ptr, am_ptr, ac_ptr, ad_ptr, rhs_d, coefx, dt)

        implicit none

        double precision, device, intent(inout) :: ap_ptr(:, :, :), am_ptr(:, :, :)
        double precision, device, intent(inout) :: ac_ptr(:, :, :), ad_ptr(:, :, :)
        double precision, device, intent(in)    :: rhs_d(:, :, :)            ! r.h.s. array
        double precision, value, intent(in)     :: coefx, dt

        integer :: i,j,k
        double precision :: inv_dmx_sub_i, inv_dmx_sub_ip

        i = (blockidx%x-1) * blockdim%x + threadidx%x
        j = (blockidx%y-1) * blockdim%y + threadidx%y
        k = (blockidx%z-1) * blockdim%z + threadidx%z

        inv_dmx_sub_i  = 1.0d0 / dmx_sub_d(i  )
        inv_dmx_sub_ip = 1.0d0 / dmx_sub_d(i+1)

        ap_ptr(j,k,i) = -coefx*inv_dmx_sub_ip*dt
        am_ptr(j,k,i) = -coefx*inv_dmx_sub_i*dt
        ac_ptr(j,k,i) =  coefx*( inv_dmx_sub_ip + inv_dmx_sub_i )*dt + 1.0d0
        ad_ptr(j,k,i) =  rhs_d(i,j,k)

    end subroutine build_LHSx_cuda

    !>
    !> @brief   Copy the solution pointer to RHS term
    !> @param   ad_ptr      Pointer to obtained solution with the linear system (device)
    !> @param   rhs_d       RHS term (device)
    !>    
    attributes(global) subroutine copy_ijk2ijk(ad_ptr, rhs_d)

        implicit none

        double precision, device, intent(in )   :: ad_ptr(:, :, :)
        double precision, device, intent(out)   ::  rhs_d(:, :, :)            ! r.h.s. array
        integer :: i,j,k

        i = (blockidx%x-1) * blockdim%x + threadidx%x
        j = (blockidx%y-1) * blockdim%y + threadidx%y
        k = (blockidx%z-1) * blockdim%z + threadidx%z

        rhs_d(i,j,k)=ad_ptr(i,j,k)

    end subroutine copy_ijk2ijk

    !>
    !> @brief   Transpose the solution pointer to RHS term
    !> @param   ad_ptr      Pointer to obtained solution with the linear system (device)
    !> @param   rhs_d       RHS term (device)
    !>    
    attributes(global) subroutine transpose_ikj2ijk(ad_ptr, rhs_d)

        implicit none

        double precision, device, intent(in )   :: ad_ptr(:, :, :)
        double precision, device, intent(out)   :: rhs_d (:, :, :)            ! r.h.s. array
        integer :: i,j,k

        i = (blockidx%x-1) * blockdim%x + threadidx%x
        j = (blockidx%y-1) * blockdim%y + threadidx%y
        k = (blockidx%z-1) * blockdim%z + threadidx%z

        rhs_d(i,j,k)=ad_ptr(i,k,j)

    end subroutine transpose_ikj2ijk

    !>
    !> @brief   Update solution
    !> @param   ad_ptr      Pointer to obtained solution with the linear system (device)
    !> @param   rhs_d       Updated solution (device)
    !>    
    attributes(global) subroutine update_theta_cuda(ad_ptr, theta_d)

        implicit none

        double precision, device, intent(in   ) :: ad_ptr(:, :, :)
        double precision, device, intent(inout) :: theta_d(0:, 0:, 0:)            ! r.h.s. array
        integer :: i,j,k

        i = (blockidx%x-1) * blockdim%x + threadidx%x
        j = (blockidx%y-1) * blockdim%y + threadidx%y
        k = (blockidx%z-1) * blockdim%z + threadidx%z

        theta_d(i,j,k) = theta_d(i,j,k) + ad_ptr(j,k,i)

    end subroutine update_theta_cuda

#endif

end module solve_theta