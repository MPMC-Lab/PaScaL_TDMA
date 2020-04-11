!======================================================================================================================
!> @file        solve_theta_plan_many.f90
!> @brief       This file contains a solver subroutine for the example problem of PaScaL_TDMA.
!> @details     The target example problem is the three-dimensional time-dependent heat conduction problem 
!>              in a unit cube domain applied with the boundary conditions of vertically constant temperature 
!>              and horizontally periodic boundaries.
!> @author      
!>              - Kiha Kim (k-kiha@yonsei.ac.kr), Department of Computational Science & Engineering, Yonsei University
!>              - Ji-Hoon Kang (jhkang@kisti.re.kr), Korea Institute of Science and Technology Information
!>              - Jung-Il Choi (jic@yonsei.ac.kr), Department of Computational Science & Engineering, Yonsei University
!>
!> @date        June 2019
!> @version     1.0
!> @par         Copyright
!>              Copyright (c) 2019 Kiha Kim and Jung-Il choi, Yonsei University and 
!>              Ji-Hoon Kang, Korea Institute of Science and Technology Information, All rights reserved.
!> @par         License     
!>              This project is release under the terms of the MIT License (see LICENSE in )
!======================================================================================================================

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

end subroutine solve_theta_plan_many

