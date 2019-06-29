!======================================================================================================================
!> @file        pascal_tdma.f90
!> @brief       PaScaL_TDMA - Parallel and Scalable Library for Tri-Diagonal Matrix Algorithm
!>              PaScal_TDMA solves many tridiagonal systems of equations in a parallel manner
!>              using a noble all-to-all communication scheme with the modified Thomas algorithm
!>              by Laszlo et al(2016).
!> @details     This module is for both of a single and many tridiagonal systems of equations.
!>              The main algorithm consists of the following five steps: 
!>              (1) Transform the original partitioned tridiagonal systems of equations 
!>                  The original partitioned tridiagonal system of equations is transformed 
!>                  to the modified tridiagonal system using the modified Thomas algorithm.
!>              (2) Build the reduced tridiagonal systems of equations
!>                  The first and last rows of the modified tridiagonal systems are assembled 
!>                  by executing the proposed all-to-all communication scheme and 
!>                  the reduced tridiagonal systems are built.
!>              (3) Solve the reduced tridiagonal systems
!>                  The reduced tridiagonal systems are solved by applying Thomas algorithm.
!>              (4) Distribute the solutions of the reduced tridiagonal systems
!>                  The solutions in Step 3 are distributed to each computing cores by executing
!>                  the inverse of the all-to-all communication in Step 2.
!>              (5) Update the solutions in the modified tridiagonal systems
!>                  The remaining unknows are updated using the modified tridiagonal systems
!>                  with the solutions obtained in Step 3 and Step 4.
!>                  
!>              Step 1 and Step 5 are similar to the method proposed by Laszlo, Gilles and Appleyard(2016)
!>              which used the parallel cyclic reduction (PCR) to build and solve the reduced tridiagonal systems.
!>              Instead of using the PCR, we develop an all-to-all communication scheme using a MPI_Ialltoall
!>              function after the modified Thomas algorithm is execued. The number of coefficients for
!>              the reduced tridiagonal systems are greatly reduced, so we can avoid the communication 
!>              bandwidth problem which is a main bottle-neck of all-to-all communications.
!>              Our algorithm is also distinguished from the work of Mattor, Williams, Hewette (1995) which
!>              assembles the undetermined coefficients of the temporary solutions in a single processor 
!>              using MPI_Gather, where load imbalances are serious.
!> 
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
!> @brief       Module for PaScaL-TDMA library 
!> @details     It contains plans for triagonal systems of equations and subroutines for solving them 
!>              using the defined plans. The operation of the library consists of the following three phases.
!>              (1) Create a data structure called a plan which has information for communication and reduced systems.
!>              (2) Solve the tridiagonal systems of equations executing from Step 1 to Step 5
!>              (3) Destroy the created plan
!>
module PaScaL_TDMA

    use mpi

    implicit none

    !> @brief   Execution plan for many tridiagonal systems of equations.
    !> @details It uses MPI_Ialltoallw function to distribute the modified tridiagonal systems to MPI processes
    !>          and build the reduced tridiagonal systems of equations. Derived datatypes are defined and  used 
    !>          for elimination of data packing and unpacking cost.
    type, public :: ptdma_plan_many

        integer :: ptdma_world          !< Single dimensional subcommunicator to assemble data for the reduced TDMA
        integer :: n_sys_rt             !< Number of tridiagonal systems to solve in each process after transpose
        integer :: n_row_rt             !< Number of rows of a reduced tridiagonal systems after transpose

        !> @{ Send buffer related variables for MPI_Ialltoallw
        integer, allocatable, dimension(:) :: ddtype_FS, count_send, displ_send
        !> @}

        !> @{ Recv. buffer related variables MPI_Ialltoallw 
        integer, allocatable, dimension(:) :: ddtype_BS, count_recv, displ_recv
        !> @}

        !> @{ Coefficient arrays after reduction, a: lower, b: diagonal, c: upper, d: rhs.
        !>    The orginal dimension (m:n) is reduced to (m:2)
        double precision, allocatable, dimension(:,:) :: A_rd, B_rd, C_rd, D_rd     
        !> @}

        !> @{ Coefficient arrays after transpose of reduced systems, a: lower, b: diagonal, c: upper, d: rhs
        !>    The reduced dimension (m:2) changes to (m/np: 2*np) after transpose.
        double precision, allocatable, dimension(:,:) :: A_rt, B_rt, C_rt, D_rt 
        !> @}

    end type ptdma_plan_many

    !> @brief   Execution plan for a single tridiagonal system of equations.
    !> @details It uses MPI_Igather function to build the reduced tridiagonal system of equations
    !>          to a specified MPI processes.
    type, public :: ptdma_plan_single

        integer :: ptdma_world          !< Single dimensional subcommunicator to assemble data for the reduced TDMA
        integer :: n_row_rt             !< Number of rows of a reduced tridiagonal system after MPI_Gather

        integer :: gather_rank          !< Destination rank of MPI_Igather
        integer :: myrank               !< Current rank ID in the communicator of ptdma_world

        !> @{ Coefficient arrays after reduction, a: lower, b: diagonal, c: upper, d: rhs.
        !>    The orginal dimension (n) is reduced to (2)
        double precision, allocatable, dimension(:) :: A_rd, B_rd, C_rd, D_rd
        !> @}

        !> @{ Coefficient arrays after transpose of a reduced system, a: lower, b: diagonal, c: upper, d: rhs
        !>    The reduced dimension (2) changes to (2*np) after transpose.
        double precision, allocatable, dimension(:) :: A_rt, B_rt, C_rt, D_rt   !< Coefficient arrays, a: lower, b: diagonal, c: upper, d: rhs
        !> @}

    end type ptdma_plan_single

    private

    public  :: PaScaL_TDMA_plan_single_create
    public  :: PaScaL_TDMA_plan_single_destroy
    public  :: PaScaL_TDMA_single_solve
    public  :: PaScaL_TDMA_single_solve_cycle

    public  :: PaScaL_TDMA_plan_many_create
    public  :: PaScaL_TDMA_plan_many_destroy
    public  :: PaScaL_TDMA_many_solve
    public  :: PaScaL_TDMA_many_solve_cycle

    contains

    !>
    !> @brief   Create a plan for a single tridiagonal system of equations
    !> @param   plan        Plan for a single tridiagonal system of equations
    !> @param   myrank      Rank ID in mpi_world
    !> @param   nprocs      Number of MPI process in mpi_world
    !> @param   mpi_world   Communicator for MPI_Gather and MPI_Scatter of a reduced system
    !> @param   gather_rank Target rank where all coeffcients are gathered to
    !>
    subroutine PaScaL_TDMA_plan_single_create(plan, myrank, nprocs, mpi_world, gather_rank)

        implicit none

        type(ptdma_plan_single), intent(inout)  :: plan
        integer, intent(in)     :: myrank, nprocs, mpi_world, gather_rank

        integer                 :: nr_rd    ! Number of rows of a reduced tridiagonal system per process, 2
        integer                 :: nr_rt    ! Number of rows of a reduced tridiagonal system after MPI_Gather

        nr_rd = 2
        nr_rt = nr_rd*nprocs

        plan%myrank = myrank
        plan%gather_rank = gather_rank
        plan%ptdma_world = mpi_world
        plan%n_row_rt = nr_rt

        allocate( plan%A_rd(1:nr_rd), plan%B_rd(1:nr_rd), plan%C_rd(1:nr_rd), plan%D_rd(1:nr_rd) )
        allocate( plan%A_rt(1:nr_rt), plan%B_rt(1:nr_rt), plan%C_rt(1:nr_rt), plan%D_rt(1:nr_rt) )
        
    end subroutine PaScaL_TDMA_plan_single_create

    !>
    !> @brief   Deallocate the allocated arrays in the defined plan_single 
    !> @param   plan        Plan for a single tridiagonal system of equations
    !>
    subroutine PaScaL_TDMA_plan_single_destroy(plan)

        implicit none

        type(ptdma_plan_single), intent(inout)  :: plan

        deallocate(plan%A_rd, plan%B_rd, plan%C_rd, plan%D_rd)
        deallocate(plan%A_rt, plan%B_rt, plan%C_rt, plan%D_rt)

    end subroutine PaScaL_TDMA_plan_single_destroy

    !>
    !> @brief   Create a plan for many tridiagonal systems of equations
    !> @param   plan        Plan for a single tridiagonal system of equations
    !> @param   n_sys       Number of tridiagonal systems of equations for process
    !> @param   myrank      Rank ID in mpi_world
    !> @param   nprocs      Number of MPI process in mpi_world
    !> @param   mpi_world   Communicator for MPI_Gather and MPI_Scatter of reduced equations
    !>
    subroutine PaScaL_TDMA_plan_many_create(plan, n_sys, myrank, nprocs, mpi_world)

        implicit none

        type(ptdma_plan_many), intent(inout)  :: plan
        integer, intent(in)     :: n_sys
        integer, intent(in)     :: myrank, nprocs, mpi_world

        integer :: i, ierr
        integer :: ista, iend                               ! First and last indices of assigned range in many tridiagonal systems of equations 
        integer :: bigsize(2), subsize(2), start(2)         ! Temporary variables of derived data type
        integer :: ns_rd, nr_rd                             ! Dimensions of many reduced tridiagonal systems
        integer :: ns_rt, nr_rt                             ! Dimensions of many reduced tridiagonal systems after transpose
        integer, allocatable, dimension(:):: ns_rt_array    ! Array specifying the number of tridiagonal systems for each process after transpose

        ! Specifying dimensions for reduced systems
        ns_rd = n_sys
        nr_rd = 2

        ! Specifying dimensions for reduced systems after transpose
        ! ns_rt         : divide the number of tridiagonal systems of equations per each process  
        ! ns_rt_array   : save the ns_rt in ns_rt_array for defining the derived data type
        ! nr_rt         : dimension of reduced tridiagonal systems in solving direction, nr_rd*nprocs
        call para_range(1, ns_rd, nprocs, myrank, ista, iend)
        ns_rt = iend - ista + 1
        allocate(ns_rt_array(0:nprocs-1))
        call MPI_Allgather(ns_rt, 1, mpi_integer, ns_rt_array,1, mpi_integer, mpi_world, ierr)
        nr_rt = nr_rd*nprocs

        ! Assign plan variables and allocate coefficient arrays
        plan%n_sys_rt = ns_rt
        plan%n_row_rt = nr_rt
        plan%ptdma_world = mpi_world

        allocate( plan%A_rd(1:ns_rd, 1:nr_rd) )
        allocate( plan%B_rd(1:ns_rd, 1:nr_rd) )
        allocate( plan%C_rd(1:ns_rd, 1:nr_rd) )
        allocate( plan%D_rd(1:ns_rd, 1:nr_rd) )
        allocate( plan%A_rt(1:ns_rt, 1:nr_rt) )
        allocate( plan%B_rt(1:ns_rt, 1:nr_rt) )
        allocate( plan%C_rt(1:ns_rt, 1:nr_rt) )
        allocate( plan%D_rt(1:ns_rt, 1:nr_rt) )

        ! Building derived data type (DDT)
        allocate(plan%ddtype_Fs(0:nprocs-1),  plan%ddtype_Bs(0:nprocs-1))

        do i=0,nprocs-1
            ! DDT for sending coefficients of the reduced tri-diagonal systems using MPI_Ialltoallw communication
            bigsize(1)=ns_rd
            bigsize(2)=nr_rd
            subsize(1)=ns_rt_array(i)
            subsize(2)=nr_rd
            start(1)=sum(ns_rt_array(0:i)) - ns_rt_array(i)
            start(2)=0
            call MPI_Type_create_subarray(  2, bigsize, subsize, start,                     &
                                            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,        &
                                            plan%ddtype_Fs(i), ierr )
            call MPI_Type_commit(plan%ddtype_Fs(i), ierr)
            ! DDT for receiving coefficients for the transposed systems of reduction using MPI_Ialltoallw communication
            bigsize(1)=ns_rt
            bigsize(2)=nr_rt
            subsize(1)=ns_rt
            subsize(2)=nr_rd
            start(1)=0
            start(2)=nr_rd*i
            call MPI_Type_create_subarray(  2, bigsize, subsize, start,                     &
                                            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,        &
                                            plan%ddtype_Bs(i), ierr )
            call MPI_Type_commit(plan%ddtype_Bs(i), ierr)
        enddo

        ! Buffer counts and displacements for MPI_Ialltoallw
        ! All Buffer counts are 1 and displacements are 0 due to the defined derived data type
        allocate(plan%count_send(0:nprocs-1), plan%displ_send(0:nprocs-1))
        allocate(plan%count_recv(0:nprocs-1), plan%displ_recv(0:nprocs-1))
        plan%count_send=1; plan%displ_send=0
        plan%count_recv=1; plan%displ_recv=0

        ! Deallocate local array
        if(allocated(ns_rt_array)) deallocate(ns_rt_array)

    end subroutine PaScaL_TDMA_plan_many_create

    !>
    !> @brief   Destroy the allocated arrays in the defined plan_many
    !> @param   plan        Plan for many tridiagonal systems of equations
    !>
    subroutine PaScaL_TDMA_plan_many_destroy(plan)

        implicit none

        type(ptdma_plan_many), intent(inout)  :: plan

        deallocate(plan%ddtype_Fs,  plan%ddtype_Bs)
        deallocate(plan%count_send, plan%displ_send)
        deallocate(plan%count_recv, plan%displ_recv)
        deallocate(plan%A_rd, plan%B_rd, plan%C_rd, plan%D_rd)
        deallocate(plan%A_rt, plan%B_rt, plan%C_rt, plan%D_rt)

    end subroutine PaScaL_TDMA_plan_many_destroy

    !>
    !> @brief   Solve a single tridiagonal system of equation
    !> @param   plan        Plan for a single tridiagonal system of equation
    !> @param   A           Coefficients in lower diagonal elements
    !> @param   B           Coefficients in diagonal elements
    !> @param   C           Coefficients in upper diagonal elements
    !> @param   D           Coefficients in right-hand side terms
    !> @param   n_row       Number of rows in each process, dimension of a tridiagonal matrix N divided by nprocs
    !>
    subroutine PaScaL_TDMA_single_solve(plan, A, B, C, D, n_row)
    
        implicit none

        type(ptdma_plan_single), intent(inout)   :: plan
        double precision,   intent(inout)   :: A(1:n_row), B(1:n_row), C(1:n_row), D(1:n_row)
        integer, intent(in) :: n_row

        ! Temporary variables for computation and parameters for MPI functions
        double precision :: r
        integer :: i
        integer :: request(4), ierr

        ! Reduction step : elimination of lower diagonal elements
        A(1) = A(1)/B(1)
        D(1) = D(1)/B(1)
        C(1) = C(1)/B(1)

        A(2) = A(2)/B(2)
        D(2) = D(2)/B(2)
        C(2) = C(2)/B(2)
    
        do i=3,n_row
            r    =  1.d0/(B(i)-A(i)*C(i-1))
            D(i) =  r*(D(i)-A(i)*D(i-1))
            C(i) =  r*C(i)
            A(i) = -r*A(i)*A(i-1)
        enddo
    
        ! Reduction step : elimination of upper diagonal elements
        do i=n_row-2,2,-1
            D(i) = D(i)-C(i)*D(i+1)
            A(i) = A(i)-C(i)*A(i+1)
            C(i) =-C(i)*C(i+1)
        enddo

        r = 1.d0/(1.d0-A(2)*C(1))
        D(1) =  r*(D(1)-C(1)*D(2))
        A(1) =  r*A(1)
        C(1) = -r*C(1)*C(2)

        ! Construct a reduced tridiagonal system of equations per each rank. Each process has two reduced rows.
        plan%A_rd(1) = A(1); plan%A_rd(2) = A(n_row)
        plan%B_rd(1) = 1.d0; plan%B_rd(2) = 1.d0
        plan%C_rd(1) = C(1); plan%C_rd(2) = C(n_row)
        plan%D_rd(1) = D(1); plan%D_rd(2) = D(n_row)

        ! Gather the coefficients of the reduced tridiagonal system to a defined rank, plan%gather_rank
        call MPI_Igather(plan%A_rd, 2, MPI_DOUBLE_PRECISION, &
                         plan%A_rt, 2, MPI_DOUBLE_PRECISION, &
                         plan%gather_rank, plan%ptdma_world, request(1), ierr)
        call MPI_Igather(plan%B_rd, 2, MPI_DOUBLE_PRECISION, &
                         plan%B_rt, 2, MPI_DOUBLE_PRECISION, &
                         plan%gather_rank, plan%ptdma_world, request(2), ierr)
        call MPI_Igather(plan%C_rd, 2, MPI_DOUBLE_PRECISION, &
                         plan%C_rt, 2, MPI_DOUBLE_PRECISION, &
                         plan%gather_rank, plan%ptdma_world, request(3), ierr)
        call MPI_Igather(plan%D_rd, 2, MPI_DOUBLE_PRECISION, &
                         plan%D_rt, 2, MPI_DOUBLE_PRECISION, &
                         plan%gather_rank, plan%ptdma_world, request(4), ierr)
        call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)

        ! Solve the reduced tridiagonal system on plan%gather_rank
        if(plan%myrank == plan%gather_rank) then
            call tdma_single(plan%A_rt,plan%B_rt,plan%C_rt,plan%D_rt, plan%n_row_rt)
        endif

        ! Scatter the solutions to each rank
        call MPI_Iscatter(plan%D_rt, 2, MPI_DOUBLE_PRECISION, &
                          plan%D_rd, 2, MPI_DOUBLE_PRECISION, &
                          plan%gather_rank, plan%ptdma_world, request(1), ierr)

        call MPI_Waitall(1, request, MPI_STATUSES_IGNORE, ierr)

        ! Update solutions of the modified tridiagonal system with the solutions of the reduced tridiagonal system.
        D(1 ) = plan%D_rd(1)
        D(n_row) = plan%D_rd(2)
        do i=2,n_row-1
            D(i) = D(i)-A(i)*D(1)-C(i)*D(n_row)
        enddo

    end subroutine PaScaL_TDMA_single_solve

    !>
    !> @brief   Solve a single cyclic tridiagonal system of equations
    !> @param   plan        Plan for a single tridiagonal system of equations
    !> @param   A           Coefficients in lower diagonal elements
    !> @param   B           Coefficients in diagonal elements
    !> @param   C           Coefficients in upper diagonal elements
    !> @param   D           Coefficients in right-hand side terms
    !> @param   n_row       Number of rows in each process, dimension of a tridiagonal matrix N divided by nprocs
    !>
    subroutine PaScaL_TDMA_single_solve_cycle(plan, A, B, C, D, n_row)

        implicit none

        type(ptdma_plan_single), intent(inout)   :: plan
        double precision, intent(inout)     :: A(1:n_row), B(1:n_row), C(1:n_row), D(1:n_row)
        integer, intent(in)                 :: n_row

        ! Temporary variables for computation and parameters for MPI functions
        integer :: i, request(4), ierr
        double precision :: rr

        ! The modified Thomas algorithm : elimination of lower diagonal elements
        A(1) = A(1)/B(1)
        D(1) = D(1)/B(1)
        C(1) = C(1)/B(1)

        A(2) = A(2)/B(2)
        D(2) = D(2)/B(2)
        C(2) = C(2)/B(2)

        do i=3,n_row
            rr = 1.d0/(B(i)-A(i)*C(i-1))
            D(i) =  rr*(D(i)-A(i)*D(i-1))
            C(i) =  rr*C(i)
            A(i) = -rr*A(i)*A(i-1)
        enddo

        ! The modified Thomas algorithm : elimination of upper diagonal elements
        do i=n_row-2,2,-1
            D(i) =  D(i)-C(i)*D(i+1)
            A(i) =  A(i)-C(i)*A(i+1)
            C(i) = -C(i)*C(i+1)
        enddo
        
        rr = 1.d0/(1.d0-A(2)*C(1))
        D(1) =  rr*(D(1)-C(1)*D(2))
        A(1) =  rr*A(1)
        C(1) = -rr*C(1)*C(2)

        ! Construct a reduced tridiagonal system of equations per each rank. Each process has two reduced rows.
        plan%A_rd(1) = A(1); plan%A_rd(2) = A(n_row)
        plan%B_rd(1) = 1.d0; plan%B_rd(2) = 1.d0
        plan%C_rd(1) = C(1); plan%C_rd(2) = C(n_row)
        plan%D_rd(1) = D(1); plan%D_rd(2) = D(n_row)

        ! Gather the coefficients of the reduced tridiagonal system to a defined rank, plan%gather_rank
        call MPI_Igather(plan%A_rd, 2, MPI_DOUBLE_PRECISION, &
                         plan%A_rt, 2, MPI_DOUBLE_PRECISION, &
                         plan%gather_rank, plan%ptdma_world, request(1), ierr)
        call MPI_Igather(plan%B_rd, 2, MPI_DOUBLE_PRECISION, &
                         plan%B_rt, 2, MPI_DOUBLE_PRECISION, &
                         plan%gather_rank, plan%ptdma_world, request(2), ierr)
        call MPI_Igather(plan%C_rd, 2, MPI_DOUBLE_PRECISION, &
                         plan%C_rt, 2, MPI_DOUBLE_PRECISION, &
                         plan%gather_rank, plan%ptdma_world, request(3), ierr)
        call MPI_Igather(plan%D_rd, 2, MPI_DOUBLE_PRECISION, &
                         plan%D_rt, 2, MPI_DOUBLE_PRECISION, &
                         plan%gather_rank, plan%ptdma_world, request(4), ierr)

        call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)

        ! Solve the reduced cyclic tridiagonal system on plan%gather_rank
        if(plan%myrank == plan%gather_rank) then
            call tdma_cycl_single(plan%A_rt,plan%B_rt,plan%C_rt,plan%D_rt, plan%n_row_rt)
        endif

        ! Distribute the solutions to each rank
        call MPI_Iscatter(plan%D_rt, 2, MPI_DOUBLE_PRECISION, &
                          plan%D_rd, 2, MPI_DOUBLE_PRECISION, &
                          plan%gather_rank, plan%ptdma_world, request(1), ierr)

        call MPI_Waitall(1, request, MPI_STATUSES_IGNORE, ierr)

        ! Update solutions of the modified tridiagonal system with the solutions of the reduced tridiagonal system.
        D(1 ) = plan%D_rd(1)
        D(n_row) = plan%D_rd(2)

        do i=2, n_row-1
            D(i) = D(i)-A(i)*D(1)-C(i)*D(n_row)
        enddo

    end subroutine PaScaL_TDMA_single_solve_cycle

    !>
    !> @brief   Solve many tridiagonal systems of equations
    !> @param   plan        Plan for many tridiagonal systems of equations
    !> @param   A           Coefficients in lower diagonal elements
    !> @param   B           Coefficients in diagonal elements
    !> @param   C           Coefficients in upper diagonal elements
    !> @param   D           Coefficients in right-hand side terms
    !> @param   n_sys       Number of tridiagonal systems per process
    !> @param   n_row       Number of rows in each process, dimension of a tridiagonal matrix N divided by nprocs
    !>
    subroutine PaScaL_TDMA_many_solve(plan, A, B, C, D, n_sys, n_row)

        implicit none

        type(ptdma_plan_many), intent(inout)   :: plan
        double precision, intent(inout)     :: A(1:n_sys,1:n_row), B(1:n_sys,1:n_row), C(1:n_sys,1:n_row), D(1:n_sys,1:n_row)
        integer, intent(in)                 :: n_sys, n_row

        ! Temporary variables for computation and parameters for MPI functions
        integer :: i, j
        integer :: request(4),ierr
        double precision :: r

        ! The modified Thomas algorithm : elimination of lower diagonal elements. 
        ! First index is to indicate a number of independent many tridiagonal systems to use vectorization
        ! Second index is to indicate a row number in a partitioned tridiagonal system 
        do i=1, n_sys
            A(i,1) = A(i,1)/B(i,1)
            D(i,1) = D(i,1)/B(i,1)
            C(i,1) = C(i,1)/B(i,1)

            A(i,2) = A(i,2)/B(i,2)
            D(i,2) = D(i,2)/B(i,2)
            C(i,2) = C(i,2)/B(i,2)
        enddo
    
        do j=3, n_row
            do i=1, n_sys
                r    =    1.d0/(B(i,j)-A(i,j)*C(i,j-1))
                D(i,j) =  r*(D(i,j)-A(i,j)*D(i,j-1))
                C(i,j) =  r*C(i,j)
                A(i,j) = -r*A(i,j)*A(i,j-1)
            enddo
        enddo
    
        ! The modified Thomas algorithm : elimination of upper diagonal elements
        do j=n_row-2, 2, -1
            do i=1, n_sys
                D(i,j) = D(i,j)-C(i,j)*D(i,j+1)
                A(i,j) = A(i,j)-C(i,j)*A(i,j+1)
                C(i,j) =-C(i,j)*C(i,j+1)
            enddo
        enddo

        do i=1, n_sys
            r = 1.d0/(1.d0-A(i,2)*C(i,1))
            D(i,1) =  r*(D(i,1)-C(i,1)*D(i,2))
            A(i,1) =  r*A(i,1)
            C(i,1) = -r*C(i,1)*C(i,2)

            ! Construct many reduced tridiagonal systems per each rank. Each process has two rows of reduced systems.
            plan%A_rd(i,1) = A(i,1); plan%A_rd(i,2) = A(i,n_row)
            plan%B_rd(i,1) = 1.d0  ; plan%B_rd(i,2) = 1.d0
            plan%C_rd(i,1) = C(i,1); plan%C_rd(i,2) = C(i,n_row)
            plan%D_rd(i,1) = D(i,1); plan%D_rd(i,2) = D(i,n_row)
        enddo

        
        ! Transpose the reduced systems of equations for TDMA using MPI_Ialltoallw and derived datatypes
        call MPI_Ialltoallw(plan%A_rd, plan%count_send, plan%displ_send, plan%ddtype_Fs, &
                            plan%A_rt, plan%count_recv, plan%displ_recv, plan%ddtype_Bs, &
                            plan%ptdma_world, request(1), ierr)
        call MPI_Ialltoallw(plan%B_rd, plan%count_send, plan%displ_send, plan%ddtype_Fs, &
                            plan%B_rt, plan%count_recv, plan%displ_recv, plan%ddtype_Bs, &
                            plan%ptdma_world, request(2), ierr)
        call MPI_Ialltoallw(plan%C_rd, plan%count_send, plan%displ_send, plan%ddtype_Fs, &
                            plan%C_rt, plan%count_recv, plan%displ_recv, plan%ddtype_Bs, &
                            plan%ptdma_world, request(3), ierr)
        call MPI_Ialltoallw(plan%D_rd, plan%count_send, plan%displ_send, plan%ddtype_Fs, &
                            plan%D_rt, plan%count_recv, plan%displ_recv, plan%ddtype_Bs, &
                            plan%ptdma_world, request(4), ierr)

        call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)

        ! Solve the reduced tridiagonal systems of equations using Thomas algorithm
        call tdma_many(plan%A_rt,plan%B_rt,plan%C_rt,plan%D_rt, plan%n_sys_rt, plan%n_row_rt)

        ! Transpose the obtained solutions to original reduced forms using MPI_Ialltoallw and derived datatypes
        call MPI_Ialltoallw(plan%D_rt, plan%count_recv, plan%displ_recv, plan%ddtype_Bs, &
                            plan%D_rd, plan%count_send, plan%displ_send, plan%ddtype_Fs, &
                            plan%ptdma_world, request(1), ierr)
        call MPI_Waitall(1, request, MPI_STATUSES_IGNORE, ierr)

        ! Update solutions of the modified tridiagonal system with the solutions of the reduced tridiagonal system.
        do i=1,n_sys
            D(i,1 ) = plan%D_rd(i,1)
            D(i,n_row) = plan%D_rd(i,2)
        enddo

        do j=2,n_row-1
            do i=1,n_sys
                D(i,j) = D(i,j)-A(i,j)*D(i,1)-C(i,j)*D(i,n_row)
            enddo
        enddo

    end subroutine PaScaL_TDMA_many_solve

    !>
    !> @brief   Solve many cyclic tridiagonal systems of equations
    !> @param   plan        Plan for many tridiagonal systems of equations
    !> @param   A           Coefficients in lower diagonal elements
    !> @param   B           Coefficients in diagonal elements
    !> @param   C           Coefficients in upper diagonal elements
    !> @param   D           Coefficients in right-hand side terms
    !> @param   n_sys       Number of tridiagonal systems per process
    !> @param   n_row       Number of rows in each process, dimension of a tridiagonal matrix N divided by nprocs
    !>
    subroutine PaScaL_TDMA_many_solve_cycle(plan, A, B, C, D, n_sys, n_row)

        implicit none

        type(ptdma_plan_many), intent(inout)   :: plan
        double precision, intent(inout)     :: A(1:n_sys,1:n_row), B(1:n_sys,1:n_row), C(1:n_sys,1:n_row), D(1:n_sys,1:n_row)
        integer, intent(in)                 :: n_sys, n_row

        ! Temporary variables for computation and parameters for MPI functions
        integer :: i,j
        integer :: request(4), ierr
        double precision :: r

        ! The modified Thomas algorithm : elimination of lower diagonal elements. 
        ! First index is to indicate a number of independent many tridiagonal systems to use vectorization
        ! Second index is to indicate a row number in a partitioned tridiagonal system 
        do i=1,n_sys
            A(i,1) = A(i,1)/B(i,1)
            D(i,1) = D(i,1)/B(i,1)
            C(i,1) = C(i,1)/B(i,1)

            A(i,2) = A(i,2)/B(i,2)
            D(i,2) = D(i,2)/B(i,2)
            C(i,2) = C(i,2)/B(i,2)
        enddo
    
        do j=3,n_row
            do i=1,n_sys
                r =    1.d0/(B(i,j)-A(i,j)*C(i,j-1))
                D(i,j) =  r*(D(i,j)-A(i,j)*D(i,j-1))
                C(i,j) =  r*C(i,j)
                A(i,j) = -r*A(i,j)*A(i,j-1)
            enddo
        enddo
    
        ! The modified Thomas algorithm : elimination of upper diagonal elements
        do j=n_row-2,2,-1
            do i=1,n_sys
                D(i,j) = D(i,j)-C(i,j)*D(i,j+1)
                A(i,j) = A(i,j)-C(i,j)*A(i,j+1)
                C(i,j) =-C(i,j)*C(i,j+1)
            enddo
        enddo

        do i=1,n_sys
            r = 1.d0/(1.d0-A(i,2)*C(i,1))
            D(i,1) =  r*(D(i,1)-C(i,1)*D(i,2))
            A(i,1) =  r*A(i,1)
            C(i,1) = -r*C(i,1)*C(i,2)

            ! Construct the reduced tridiagonal equations per each rank. Each process has two rows of reduced systems.
            plan%A_rd(i,1) = A(i,1); plan%A_rd(i,2) = A(i,n_row)
            plan%B_rd(i,1) = 1.d0  ; plan%B_rd(i,2) = 1.d0
            plan%C_rd(i,1) = C(i,1); plan%C_rd(i,2) = C(i,n_row)
            plan%D_rd(i,1) = D(i,1); plan%D_rd(i,2) = D(i,n_row)
        enddo
        
        ! Transpose the reduced systems of equations for TDMA using MPI_Ialltoallw and derived datatypes
        call MPI_Ialltoallw(plan%A_rd, plan%count_send, plan%displ_send, plan%ddtype_Fs, &
                            plan%A_rt, plan%count_recv, plan%displ_recv, plan%ddtype_Bs, &
                            plan%ptdma_world, request(1), ierr)
        call MPI_Ialltoallw(plan%B_rd, plan%count_send, plan%displ_send, plan%ddtype_Fs, &
                            plan%B_rt, plan%count_recv, plan%displ_recv, plan%ddtype_Bs, &
                            plan%ptdma_world, request(2), ierr)
        call MPI_Ialltoallw(plan%C_rd, plan%count_send, plan%displ_send, plan%ddtype_Fs, &
                            plan%C_rt, plan%count_recv, plan%displ_recv, plan%ddtype_Bs, &
                            plan%ptdma_world, request(3), ierr)
        call MPI_Ialltoallw(plan%D_rd, plan%count_send, plan%displ_send, plan%ddtype_Fs, &
                            plan%D_rt, plan%count_recv, plan%displ_recv, plan%ddtype_Bs, &
                            plan%ptdma_world, request(4), ierr)

        call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)

        ! Solve the reduced cyclic tridiagonal systems of equations using cyclic TDMA
        call tdma_cycl_many(plan%A_rt, plan%B_rt, plan%C_rt, plan%D_rt, plan%n_sys_rt, plan%n_row_rt)

        ! Transpose the obtained solutions to original reduced forms using MPI_Ialltoallw and derived datatypes
        call MPI_Ialltoallw(plan%D_rt, plan%count_recv, plan%displ_recv, plan%ddtype_Bs, &
                            plan%D_rd, plan%count_send, plan%displ_send, plan%ddtype_Fs, &
                            plan%ptdma_world, request(1), ierr)
        call MPI_Waitall(1, request, MPI_STATUSES_IGNORE, ierr)

        ! Update solutions of the modified tridiagonal system with the solutions of the reduced tridiagonal system.
        do i=1,n_sys
            D(i,1 ) = plan%D_rd(i,1)
            D(i,n_row) = plan%D_rd(i,2)
        enddo

        do j=2,n_row-1
            do i=1,n_sys
                D(i,j) = D(i,j)-A(i,j)*D(i,1)-C(i,j)*D(i,n_row)
            enddo
        enddo

    end subroutine PaScaL_TDMA_many_solve_cycle

end module PaScaL_TDMA