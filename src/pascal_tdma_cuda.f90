!======================================================================================================================
!> @file        pascal_tdma_cuda.cuf
!> @brief       PaScaL_TDMA - Parallel and Scalable Library for TriDiagonal Matrix Algorithm
!> @details     PaScaL_TDMA includes a CUDA implementation of PaScaL_TDMA, which accelerates 
!>              to solve many tridiagonal systems in multi-dimensional partial differential equations on GPU.
!>              It adopts the pipeline copy within the shared memory for the forward elemination and 
!>              backward substitution procudures of TDMA to reduce global memory access.
!>              For the main algorithm of PaScaL_TDMA, see also https://github.com/MPMC-Lab/PaScaL_TDMA.
!> 
!> @author      
!>              - Mingyu Yang (yang926@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
!>              - Ji-Hoon Kang (jhkang@kisti.re.kr), Korea Institute of Science and Technology Information
!>              - Ki-Ha Kim (k-kiha@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
!>              - Jung-Il Choi (jic@yonsei.ac.kr), School of Mathematics and Computing (Computational Science and Engineering), Yonsei University
!>
!> @date        May 2023
!> @version     2.0
!> @par         Copyright
!>              Copyright (c) 2019-2023 Mingyu Yang, Ki-Ha Kim and Jung-Il choi, Yonsei University and 
!>              Ji-Hoon Kang, Korea Institute of Science and Technology Information, All rights reserved.
!> @par         License     
!>              This project is release under the terms of the MIT License (see LICENSE file).
!======================================================================================================================

!>
!> @brief       Module for PaScaL_TDMA library with CUDA.
!> @details     It contains plans for tridiagonal systems of equations and subroutines for solving them 
!>              using the defined plans. The operation of the library includes the following three phases:
!>              (1) Create a data structure called a plan with the information for communication and reduced systems.
!>              (2) Solve the tridiagonal systems of equations executing from Step 1 to Step 5
!>              (3) Destroy the created plan
!>
module PaScaL_TDMA_cuda

    use mpi
    use tdma
    use cudafor

    implicit none

    !> @brief   Execution plan for many tridiagonal systems of equations.
    !> @details It uses MPI_Ialltoall function to distribute the modified tridiagonal systems to MPI processes
    !>          and build the reduced tridiagonal systems of equations. Currently it supports the equal size of domains.
	!>          It does not use derived datatypes.
    
    type, public :: ptdma_plan_many_cuda

        private

        integer :: ptdma_world      !< Single dimensional subcommunicator to assemble data for the reduced TDMA
        integer :: nprocs     		!< Communicator size of ptdma_world

        integer :: nx_sys           !< Number of tridiagonal systems of equations in x-direction per process
        integer :: ny_sys           !< Number of tridiagonal systems of equations in y-direction per process
        integer :: nz_row           !< Row size of partitioned tridiagonal matrix in z-direction per process

        integer :: nx_sys_rt        !< Number of tridiagonal systems to be solved in each process after transpose
        integer :: nz_row_rt        !< Number of rows of a reduced tridiagonal systems after transpose
        integer :: nz_row_rd        !< Number of rows of a reduced tridiagonal systems before transpose, 2.

        !> @{ Coefficient arrays after reduction, a: lower, b: diagonal, c: upper, d: rhs.
        !>    The orginal dimension (m:n) is reduced to (m:2)
        !>    The arrays are allocated in the device memory.
        double precision, allocatable, device, dimension(:,:,:) :: a_rd_d, b_rd_d, c_rd_d, d_rd_d
        !> @}

        !> @{ Coefficient arrays for cyclic TDMA in case of nprocs = 1.
        !>    The arrays are allocated in the device memory.
        double precision, allocatable, device, dimension(:,:,:) :: e_buff
        !> @}

        !> @{ Coefficient arrays after transpose of reduced systems, a: lower, b: diagonal, c: upper, d: rhs
        !>    The reduced dimension (m:2) changes to (m/np: 2*np) after transpose.
        !>    The arrays are allocated in the device memory.
        double precision, allocatable, device, dimension(:,:,:) :: a_rt_d, b_rt_d, c_rt_d, d_rt_d, e_rt_d
        !> @}

        !> @{ Buffers allocated in device memory for all-to-all communication
        double precision, allocatable, device, dimension(:) :: sendbuf, recvbuf
        !> @}

        !> @{ Threads and blocks for CUDA kernel
        type(dim3)  :: threads, blocks, blocks_rt, blocks_alltoall
        !> @

        integer     :: shared_buffer_size   !< shared buffer size
    
    end type ptdma_plan_many_cuda

    private
    public  :: PaScaL_TDMA_plan_many_create_cuda
    public  :: PaScaL_TDMA_plan_many_destroy_cuda
    public  :: PaScaL_TDMA_many_solve_cuda
    public  :: PaScaL_TDMA_many_solve_cycle_cuda 

    contains

    !>
    !> @brief   Create a plan for many tridiagonal systems of equations.
    !> @param   p           Plan for a many tridiagonal system of equations
    !> @param   nx_sys      Number of tridiagonal systems of equations in x-direction per process
    !> @param   ny_sys      Number of tridiagonal systems of equations in y-direction per process
    !> @param   nz_row      Row size of partitioned tridiagonal matrix in z-direction per process
    !> @param   myrank      Rank ID in mpi_world
    !> @param   nprocs      Number of MPI process in mpi_world
    !> @param   mpi_world   Communicator for MPI_Gather and MPI_Scatter of reduced equations
    !>
    subroutine PaScaL_TDMA_plan_many_create_cuda(p, nx_sys, ny_sys, nz_row, myrank, nprocs, mpi_world, thread_in)

        implicit none
        
        type(ptdma_plan_many_cuda), intent(inout)  :: p
        type(dim3) :: thread_in

        integer, intent(in)     :: nx_sys
        integer, intent(in)     :: ny_sys
        integer, intent(in)     :: nz_row
        integer, intent(in)     :: myrank, nprocs, mpi_world

        integer :: nx_block, ny_block, nx_rt_block
        integer :: i, ierr
        integer :: ista, iend           ! First and last row indices of many tridiagonal systems of equations 
        integer :: nx_sys_rd, nz_row_rd ! Dimensions of many reduced tridiagonal systems
        integer :: nx_sys_rt, nz_row_rt ! Dimensions of many reduced tridiagonal systems after transpose

        ! Assign plan variables and allocate coefficient arrays.
        p%nx_sys        = nx_sys
        p%ny_sys        = ny_sys
        p%nz_row        = nz_row
        p%ptdma_world   = mpi_world
        p%nprocs  = nprocs
    
        ! Specify dimensions for reduced systems.
        nx_sys_rd = nx_sys
        nz_row_rd = 2

        ! Specify dimensions for reduced systems after transpose.
        ! nx_sys_rt         : divide the number of tridiagonal systems of equations per each process  
        ! nx_sys_rt_array   : save the nx_sys_rt in nx_sys_rt_array for defining the DDT
        ! nz_row_rt         : dimensions of the reduced tridiagonal systems in the solving direction, nz_row_rd*nprocs
        call para_range(1, nx_sys_rd, nprocs, myrank, ista, iend)
        nx_sys_rt = iend - ista + 1
        nz_row_rt = nz_row_rd*nprocs
        
        ! Specify dimensions for reduced systems.
        p%nx_sys_rt     = nx_sys_rt
        p%nz_row_rt     = nz_row_rt
        p%nz_row_rd     = nz_row_rd

        ! Allocate coefficient arrays.
        if(p%nprocs.eq.1) then
            allocate( p%e_buff(nx_sys, ny_sys, nz_row) );  p%e_buff = 0.0d0
        else
            allocate( p%a_rd_d(nx_sys_rd, ny_sys, nz_row_rd) );  p%a_rd_d  = 0.0d0
            allocate( p%b_rd_d(nx_sys_rd, ny_sys, nz_row_rd) );  p%b_rd_d  = 0.0d0
            allocate( p%c_rd_d(nx_sys_rd, ny_sys, nz_row_rd) );  p%c_rd_d  = 0.0d0
            allocate( p%d_rd_d(nx_sys_rd, ny_sys, nz_row_rd) );  p%d_rd_d  = 0.0d0
            allocate( p%a_rt_d(nx_sys_rt, ny_sys, nz_row_rt) );  p%a_rt_d  = 0.0d0
            allocate( p%b_rt_d(nx_sys_rt, ny_sys, nz_row_rt) );  p%b_rt_d  = 0.0d0
            allocate( p%c_rt_d(nx_sys_rt, ny_sys, nz_row_rt) );  p%c_rt_d  = 0.0d0
            allocate( p%d_rt_d(nx_sys_rt, ny_sys, nz_row_rt) );  p%d_rt_d  = 0.0d0
            allocate( p%e_rt_d(nx_sys_rt, ny_sys, nz_row_rt) );  p%e_rt_d  = 0.0d0
            allocate( p%sendbuf(nx_sys*ny_sys*nz_row_rd))     ;  p%sendbuf = 0.0d0
            allocate( p%recvbuf(nx_sys*ny_sys*nz_row_rd))     ;  p%recvbuf = 0.0d0
        endif

        ! Check the thread and block
        nx_block = nx_sys/thread_in%x
        if(nx_block.eq.0 .or. (mod(nx_sys, thread_in%x).ne.0)) then
            print '(a,i5,a,i5)', '[Error] nx_sys should be a multiple of thread_in%x. &
                                    thread_in%x = ',thread_in%x,', nx_sys = ',nx_sys
            call MPI_Finalize(ierr)
            stop
        endif

        nx_rt_block = nx_sys_rt/thread_in%x
        if(nx_rt_block.eq.0 .or. (mod(nx_sys_rt, thread_in%x).ne.0)) then
            print '(a,i5,a,i5)', '[Error] nx_sys_rt should be a multiple of thread_in%x. &
                                    thread_in%x = ',thread_in%x,', nx_sys_rt = ', nx_sys_rt
            call MPI_Finalize(ierr)
            stop
        endif

        ny_block = ny_sys/thread_in%y
        if(ny_block.eq.0 .or. (mod(ny_sys, thread_in%y).ne.0)) then
            print '(a,i5,a,i5)', '[Error] ny_sys should be a multiple of thread_in%y. &
                                    thread_in%y = ',thread_in%y,', ny_sys = ',ny_sys
            call MPI_Finalize(ierr)
            stop
        endif

        ! Define threads and blocks of dim3 type
        p%threads         = dim3( thread_in%x, thread_in%y,    1)
        p%blocks          = dim3( nx_block,    ny_block,       1)
        p%blocks_rt       = dim3( nx_rt_block, ny_block,       1)
        p%blocks_alltoall = dim3( nx_rt_block, ny_block,       nz_row_rd)

        ! Define the buffer size of shared memory for pipeline copy
        p%shared_buffer_size = kind(0.0d0)*(1+thread_in%x)*thread_in%y

    end subroutine PaScaL_TDMA_plan_many_create_cuda

    !>
    !> @brief   Destroy the allocated arrays in the defined plan_many.
    !> @param   p           Plan for many tridiagonal systems of equations
    !>
    subroutine PaScaL_TDMA_plan_many_destroy_cuda(p)
        implicit none

        type(ptdma_plan_many_cuda), intent(inout)  :: p

        if(p%nprocs.eq.1) then
            deallocate(p%e_buff) 
        else
            deallocate(p%a_rd_d, p%b_rd_d, p%c_rd_d, p%d_rd_d)
            deallocate(p%a_rt_d, p%b_rt_d, p%c_rt_d, p%d_rt_d, p%e_rt_d)
            deallocate(p%sendbuf, p%recvbuf)
        endif       

    end subroutine PaScaL_TDMA_plan_many_destroy_cuda

    !>
    !> @brief   Solve many tridiagonal systems of equations.
    !> @param   p           Plan for many tridiagonal systems of equations
    !> @param   a_d         Coefficient array of lower diagonal elements
    !> @param   b_d         Coefficient array of diagonal elements
    !> @param   c_d         Coefficient array of upper diagonal elements
    !> @param   d_d         Coefficient array of right-hand side terms
    !>
    subroutine PaScaL_TDMA_many_solve_cuda(p, a_d, b_d, c_d, d_d)

        implicit none

        type(ptdma_plan_many_cuda), intent(inout)   :: p
        double precision, device, intent(inout)     :: a_d(:, :, :), b_d(:, :, :)
        double precision, device, intent(inout)     :: c_d(:, :, :), d_d(:, :, :)

        if(p%nprocs.eq.1) then
            ! Solve the tridiagonal system directly when nprocs = 1. 
            ! The size of shared memory is specified using shared_buffer_size
            call tdma_many_cuda <<<p%blocks, p%threads, 6*p%shared_buffer_size>>> & ! 6 * 8 byte ?
                                (a_d, b_d, c_d, d_d, p%nz_row)
        else
            ! The modified Thomas algorithm
            ! The size of shared memory is specified using 'shared_buffer_size'
            call PaScaL_TDMA_many_modified_Thomas_cuda <<<p%blocks, p%threads, 9*p%shared_buffer_size>>> &
                    (a_d, b_d, c_d, d_d, p%a_rd_d, p%b_rd_d, p%c_rd_d, p%d_rd_d, p%nz_row)

            ! Transpose the reduced systems of equations for TDMA using MPI_alltoall and CUDA-aware-MPI.
            call transpose_slab_xy_to_yz(p, p%a_rd_d, p%a_rt_d)
            call transpose_slab_xy_to_yz(p, p%b_rd_d, p%b_rt_d)
            call transpose_slab_xy_to_yz(p, p%c_rd_d, p%c_rt_d)
            call transpose_slab_xy_to_yz(p, p%d_rd_d, p%d_rt_d)

            ! Solve the reduced tridiagonal systems of equations using Thomas algorithm.
            ! The size of shared memory is specified using 'shared_buffer_size'
            call tdma_many_cuda <<<p%blocks_rt, p%threads, 6*p%shared_buffer_size>>> &
                                (p%a_rt_d,p%b_rt_d,p%c_rt_d,p%d_rt_d, p%nz_row_rt)

            ! Transpose the obtained solutions to original reduced forms using MPI_Alltoall and CUDA-aware-MPI.
            call transpose_slab_yz_to_xy(p, p%d_rt_d, p%d_rd_d)

            ! Update solutions of the modified tridiagonal system with the solutions of the reduced tridiagonal system.
            ! The size of shared memory is specified using 'shared_buffer_size'
            call PaScaL_TDMA_many_update_solution_cuda<<<p%blocks, p%threads, 2*p%shared_buffer_size>>> &
                                    (a_d, c_d, d_d, p%d_rd_d, p%nz_row)
        endif

    end subroutine PaScaL_TDMA_many_solve_cuda

    !>
    !> @brief   Solve many cyclic tridiagonal systems of equations.
    !> @param   p           Plan for many tridiagonal systems of equations
    !> @param   a_d         Coefficient array of lower diagonal elements
    !> @param   b_d         Coefficient array of diagonal elements
    !> @param   c_d         Coefficient array of upper diagonal elements
    !> @param   d_d         Coefficient array of right-hand side terms
    !>
    subroutine PaScaL_TDMA_many_solve_cycle_cuda(p, a_d, b_d, c_d, d_d)

        implicit none

        type(ptdma_plan_many_cuda), intent(inout)   :: p
        double precision, device, intent(inout)     :: a_d(:, :, :), b_d(:, :, :)
        double precision, device, intent(inout)     :: c_d(:, :, :), d_d(:, :, :)

        if(p%nprocs.eq.1) then
            ! Solve the cyclic tridiagonal system directly when nprocs = 1. 
            ! The size of shared memory is specified using shared_buffer_size
            call tdma_cycl_many_cuda<<<p%blocks, p%threads, 8*p%shared_buffer_size>>> &
                                    (a_d, b_d, c_d, d_d, p%e_buff, p%nz_row)
        else
            ! The modified Thomas algorithm
            ! The size of shared memory is specified using 'shared_buffer_size'
            call PaScaL_TDMA_many_modified_Thomas_cuda<<<p%blocks, p%threads, 9*p%shared_buffer_size>>> &
                (a_d, b_d, c_d, d_d, p%a_rd_d, p%b_rd_d, p%c_rd_d, p%d_rd_d, p%nz_row)
 
            ! Transpose the reduced systems of equations for TDMA using MPI_alltoall and CUDA-aware-MPI.
            call transpose_slab_xy_to_yz(p, p%a_rd_d, p%a_rt_d)
            call transpose_slab_xy_to_yz(p, p%b_rd_d, p%b_rt_d)
            call transpose_slab_xy_to_yz(p, p%c_rd_d, p%c_rt_d)
            call transpose_slab_xy_to_yz(p, p%d_rd_d, p%d_rt_d)
    
            ! Solve the reduced tridiagonal systems of equations using Thomas algorithm.
            ! The size of shared memory is specified using 'shared_buffer_size'
            call tdma_cycl_many_cuda<<<p%blocks_rt, p%threads, 8*p%shared_buffer_size>>> &
                                    (p%a_rt_d,p%b_rt_d,p%c_rt_d,p%d_rt_d,p%e_rt_d, p%nz_row_rt)

            ! Transpose the obtained solutions to original reduced forms using MPI_alltoall and CUDA-aware-MPI.
            call transpose_slab_yz_to_xy(p, p%d_rt_d, p%d_rd_d)

            ! Update solutions of the modified tridiagonal system with the solutions of the reduced tridiagonal system.
            ! The size of shared memory is specified using 'shared_buffer_size'
            call PaScaL_TDMA_many_update_solution_cuda<<<p%blocks, p%threads, 2*p%shared_buffer_size>>> &
                                                    (a_d, c_d, d_d, p%d_rd_d, p%nz_row)
        endif

    end subroutine PaScaL_TDMA_many_solve_cycle_cuda

    !>
    !> @brief   The modified Thomas algorithm : elimination of lower diagonal elements
    !> @param   a           Coefficient array of lower diagonal elements
    !> @param   b           Coefficient array of diagonal elements
    !> @param   c           Coefficient array of upper diagonal elements
    !> @param   d           Coefficient array of right-hand side terms
    !> @param   a_rd        Reduced coefficient array of lower diagonal elements
    !> @param   b_rd        Reduced coefficient array of diagonal elements
    !> @param   c_rd        Reduced coefficient array of upper diagonal elements
    !> @param   d_rd        Reduced doefficient array of right-hand side terms
    !> @param   nz_row      Row size of partitioned tridiagonal matrix in z-direction per process
    !>
    attributes(global) subroutine PaScaL_TDMA_many_modified_Thomas_cuda(a, b, c, d, a_rd, b_rd, c_rd, d_rd, nz_row)

        implicit none

        integer, value, intent(in)      :: nz_row
        double precision, device, intent(inout) :: a(:, :, :), c(:, :, :), d(:, :, :)
        double precision, device, intent(in)    :: b(:, :, :)
        double precision, device, intent(inout) :: a_rd(:, :, :), b_rd(:, :, :)
        double precision, device, intent(inout) :: c_rd(:, :, :), d_rd(:, :, :)

        ! Temporary variables for computation
        integer :: i, j, k
        integer :: ti, tj, tk
        double precision :: r

        ! Block shared memory for pipeline copy
        double precision, shared :: a1(blockdim%x+1, blockdim%y), a0(blockdim%x+1, blockdim%y)
        double precision, shared :: b1(blockdim%x+1, blockdim%y), b0(blockdim%x+1, blockdim%y)
        double precision, shared :: c1(blockdim%x+1, blockdim%y), c0(blockdim%x+1, blockdim%y)
        double precision, shared :: d1(blockdim%x+1, blockdim%y), d0(blockdim%x+1, blockdim%y)
        double precision, shared :: r0(blockdim%x+1, blockdim%y)

        ! Global index
        j = (blockidx%x-1) * blockdim%x + threadidx%x
        k = (blockidx%y-1) * blockdim%y + threadidx%y

        ! Local index in block
        tj = threadidx%x
        tk = threadidx%y

        ! The modified Thomas algorithm : elimination of lower diagonal elements. 
        ! First & second indices indicate a number of independent many tridiagonal systems for parallezation.
        ! Third index indicates a row number in a partitioned tridiagonal system .
		! Therefore, first & second indices are for thread IDs.
        ! The modified Thomas algorithm : elimination of lower diagonal elements. 

        a0(tj, tk) = a(j, k, 1)
        b0(tj, tk) = b(j, k, 1)
        c0(tj, tk) = c(j, k, 1)
        d0(tj, tk) = d(j, k, 1)

        a0(tj, tk) = a0(tj, tk) / b0(tj, tk)
        c0(tj, tk) = c0(tj, tk) / b0(tj, tk)
        d0(tj, tk) = d0(tj, tk) / b0(tj, tk)

        a(j, k,1) = a0(tj, tk)
        c(j, k,1) = c0(tj, tk)
        d(j, k,1) = d0(tj, tk)

        a1(tj, tk) = a(j, k, 2)
        b1(tj, tk) = b(j, k, 2)
        c1(tj, tk) = c(j, k, 2)
        d1(tj, tk) = d(j, k, 2)

        a1(tj, tk) = a1(tj, tk) / b1(tj, tk)
        c1(tj, tk) = c1(tj, tk) / b1(tj, tk)
        d1(tj, tk) = d1(tj, tk) / b1(tj, tk)

        a(j, k,2) = a1(tj, tk)
        c(j, k,2) = c1(tj, tk)
        d(j, k,2) = d1(tj, tk)

        do i = 3, nz_row

            ! Pipeline copy of (i-1)th data using shared memory
            a0(tj, tk) = a1(tj, tk)
            c0(tj, tk) = c1(tj, tk)
            d0(tj, tk) = d1(tj, tk)

            ! Load i-th data from global memory
            a1(tj, tk) = a(j, k,i)
            b1(tj, tk) = b(j, k,i)
            c1(tj, tk) = c(j, k,i)
            d1(tj, tk) = d(j, k,i)

            r0(tj, tk) =  1.0d0 / (b1(tj, tk) - a1(tj, tk) * c0(tj, tk))
            d1(tj, tk) =  r0(tj, tk) * (d1(tj, tk) - a1(tj, tk) * d0(tj, tk))
            c1(tj, tk) =  r0(tj, tk) * c1(tj, tk)
            a1(tj, tk) = -r0(tj, tk) * a1(tj, tk) * a0(tj, tk)

            ! Save updated i-th data to global memory
            a(j, k,i) = a1(tj, tk)
            c(j, k,i) = c1(tj, tk)
            d(j, k,i) = d1(tj, tk)

        enddo

        ! Construct many reduced tridiagonal systems per each rank. Each process has two rows of reduced systems.
        a_rd(j, k,2) = a1(tj, tk)
        b_rd(j, k,2) = 1.0d0
        c_rd(j, k,2) = c1(tj, tk)
        d_rd(j, k,2) = d1(tj, tk)

        ! Pipeline copy between arrays in shared memory
        a1(tj, tk) = a0(tj, tk)
        c1(tj, tk) = c0(tj, tk)
        d1(tj, tk) = d0(tj, tk)

        ! The modified Thomas algorithm : elimination of upper diagonal elements.
        do i = nz_row-2, 2, -1

            ! Load i-th data from global memory
            a0(tj, tk) = a(j, k,i)
            c0(tj, tk) = c(j, k,i)
            d0(tj, tk) = d(j, k,i)

            d0(tj, tk) = d0(tj, tk) - c0(tj, tk)*d1(tj, tk)
            a0(tj, tk) = a0(tj, tk) - c0(tj, tk)*a1(tj, tk)
            c0(tj, tk) =-c0(tj, tk) * c1(tj, tk)

            ! Pipeline copy of updated i-th data using shared memory
            a1(tj, tk) = a0(tj, tk)
            c1(tj, tk) = c0(tj, tk)
            d1(tj, tk) = d0(tj, tk)

            ! Save updated data to global memory
            a(j, k,i) = a0(tj, tk)
            c(j, k,i) = c0(tj, tk)
            d(j, k,i) = d0(tj, tk)
        enddo

        a0(tj, tk) = a(j, k,1)
        c0(tj, tk) = c(j, k,1)
        d0(tj, tk) = d(j, k,1)

        r0(tj, tk) = 1.0d0 / (1.0d0 - a1(tj, tk) * c0(tj, tk))
        d0(tj, tk) =  r0(tj, tk) * (d0(tj, tk) - c0(tj, tk) * d1(tj, tk))
        a0(tj, tk) =  r0(tj, tk) * a0(tj, tk)
        c0(tj, tk) = -r0(tj, tk) * c0(tj, tk) * c1(tj, tk)

        d(j, k,1) = d0(tj, tk)
        a(j, k,1) = a0(tj, tk)
        c(j, k,1) = c0(tj, tk)

        ! Construct many reduced tridiagonal systems per each rank. Each process has two rows of reduced systems.
        a_rd(j, k,1) = a0(tj, tk)
        b_rd(j, k,1) = 1.0d0
        c_rd(j, k,1) = c0(tj, tk)
        d_rd(j, k,1) = d0(tj, tk)
        

    end subroutine PaScaL_TDMA_many_modified_Thomas_cuda

    !>
    !> @brief   The modified Thomas algorithm : elimination of lower diagonal elements
    !> @param   a           Coefficient array of lower diagonal elements
    !> @param   c           Coefficient array of upper diagonal elements
    !> @param   d           Coefficient array of solution terms
    !> @param   d_rd        Reduced coefficient array of solution terms
    !> @param   nz_row      Row size of partitioned tridiagonal matrix in z-direction per process
    !>
    attributes(global) subroutine PaScaL_TDMA_many_update_solution_cuda(a, c, d, d_rd, nz_row)

        implicit none

        integer, value, intent(in)      :: nz_row
        double precision, device, intent(in)    :: a(:, :, :), c(:, :, :), d_rd(:, :, :)
        double precision, device, intent(inout) :: d(:, :, :)

        ! Temporary variables for computation
        integer :: i, j, k
        integer :: tj, tk

        ! Block shared memory
        double precision, shared :: ds(blockdim%x + 1, blockdim%y), de(blockdim%x + 1, blockdim%y)

        ! Global index
        j = (blockidx%x - 1) * blockdim%x + threadidx%x
        k = (blockidx%y - 1) * blockdim%y + threadidx%y

        ! Local index in block
        tj = threadidx%x
        tk = threadidx%y

        ! Using shared memory
		! First and second indices are for thread IDs
        ds(tj, tk) = d_rd(j, k, 1)
        de(tj, tk) = d_rd(j, k, 2)
        call syncthreads()

        ! Update solutions of the modified tridiagonal system with the solutions of the reduced tridiagonal system.
        d(j, k, 1)      = ds(tj, tk)
        d(j, k, nz_row) = de(tj, tk)

        do i = 2, nz_row-1
            d(j, k, i) = d(j, k, i) - a(j, k, i) * ds(tj, tk) - c(j, k, i) * de(tj, tk)
        enddo

    end subroutine PaScaL_TDMA_many_update_solution_cuda


    !>
    !> @brief   Subroutine to transpose x-y slab to y-z slab for solving TDM in z-direction
    !> @param   slab_xy     Coefficient array in the shape of x-y slab
    !> @param   slab_yz     Coefficient array in the shape of y-z slab
    !>
    subroutine transpose_slab_xy_to_yz (p, slab_xy, slab_yz)

        implicit none

        type(ptdma_plan_many_cuda), intent(inout)   :: p
        double precision, device, intent(in )       :: slab_xy(:, :, :)
        double precision, device, intent(out)       :: slab_yz(:, :, :)
        
        integer  :: nblksize
        integer  :: ierr
  
        nblksize = p%nx_sys * p%ny_sys * p%nz_row_rd / p%nprocs

        call mem_detach_slab_xy <<<p%blocks_alltoall, p%threads>>> &
                                (slab_xy, p%sendbuf, p%nx_sys, p%ny_sys, p%nz_row_rd, p%nprocs)

        !----- alltoall communication of sbuf to rbuf
        ierr = cudaStreamSynchronize()
        call MPI_Alltoall(p%sendbuf, nblksize, MPI_DOUBLE_PRECISION, &
                          p%recvbuf, nblksize, MPI_DOUBLE_PRECISION, &
                          p%ptdma_world, ierr)

        call mem_unite_slab_yz<<<p%blocks_alltoall, p%threads>>> &
                                (p%recvbuf, slab_yz, p%nx_sys, p%ny_sys, p%nz_row_rd, p%nprocs)
  
    end subroutine transpose_slab_xy_to_yz

    !>
    !> @brief   Subroutine to rearrange x-y slab to a 1D array for MPI_Alltoall
    !> @param   slab_xy     Coefficient array in the shape of x-y slab
    !> @param   array1D     1-D Coefficient array
    !> @param   n1, n2, n3  Dimension of array 'slab_xy'
    !> @param   nprocs      Number of processes in 'mpi_comm'
    !>
    attributes(global) subroutine mem_detach_slab_xy(slab_xy, array1D, n1, n2, n3, nprocs)

        implicit none

        integer, value, intent(in)  :: n1, n2, n3, nprocs
        double precision, device, intent(in)    :: slab_xy(:, :, :)
        double precision, device, intent(out)   :: array1D(:)

        ! Variables to calculate indices for in and out arrays
        integer :: i, j, k, kblk, n1blksize, blksize
        integer :: pos, pos_i, pos_j, pos_k

        n1blksize  = n1 / nprocs
        blksize    = n1blksize * n2 * n3

        i = (blockidx%x - 1) * blockdim%x + threadidx%x
        j = (blockidx%y - 1) * blockdim%y + threadidx%y
        k = (blockidx%z - 1) * blockdim%z + threadidx%z

        pos_i = (i - 1) * n2 * n3
        pos_j = (j - 1) * n3

        do kblk = 1, nprocs
            pos_k = k + (kblk - 1) * blksize
            pos   = pos_k + pos_j + pos_i
            array1D(pos) = slab_xy(i + (kblk - 1) * n1blksize, j, k)
        enddo

    end subroutine mem_detach_slab_xy

    !>
    !> @brief   Subroutine to rearrange a 1D array to y-z slab after MPI_Alltoall
    !> @param   array1D     1-D Coefficient array
    !> @param   slab_yz     Coefficient array in the shape of y-z slab
    !> @param   n1, n2, n3  Dimension of array 'slab_xy'
    !> @param   nprocs      Number of processes in 'mpi_comm'
    !>
    attributes(global) subroutine mem_unite_slab_yz(array1D, slab_yz, n1, n2, n3, nprocs)

        implicit none

        integer, value, intent(in)  :: n1, n2, n3, nprocs
        double precision, device, intent(in)    :: array1D(:)
        double precision, device, intent(out)   :: slab_yz(:, :, :)

        ! Variables to calculate indices for in and out arrays
        integer :: i, j, k, kblk, blksize
        integer :: pos, pos_i, pos_j, pos_k

        blksize  = n1 * n2 * n3 / nprocs

        i = (blockidx%x - 1) * blockdim%x + threadidx%x
        j = (blockidx%y - 1) * blockdim%y + threadidx%y
        k = (blockidx%z - 1) * blockdim%z + threadidx%z

        pos_i = (i - 1) * n2 * n3
        pos_j = (j - 1) * n3

        do kblk = 1, nprocs
            pos_k = k + (kblk - 1) * blksize
            pos   = pos_k + pos_j + pos_i
            slab_yz(i, j, k + (kblk - 1) * n3) = array1D(pos)
        enddo

    end subroutine mem_unite_slab_yz

    !>
    !> @brief   Subroutine to transpose y-z slab (d_rt_d) to x-y slab (d_rd_d) after solving TDM in z-direction
    !> @param   p           Plan for a single tridiagonal system of equations
    !> @param   slab_yz     Coefficient array in the shape of y-z slab
    !> @param   slab_xy     Coefficient array in the shape of x-y slab
    !>
    subroutine transpose_slab_yz_to_xy(p, slab_yz, slab_xy)

        implicit none

        type(ptdma_plan_many_cuda), intent(inout)      :: p
        double precision, device, intent(in )   :: slab_yz(:, :, :)
        double precision, device, intent(out)   :: slab_xy(:, :, :)

        integer  :: nblksize
        integer  :: ierr
  
        nblksize   = p%nx_sys * p%ny_sys * p%nz_row_rd / p%nprocs

        call mem_detach_slab_yz<<<p%blocks_alltoall, p%threads>>> &
                                (slab_yz, p%sendbuf, p%nx_sys, p%ny_sys, p%nz_row_rd, p%nprocs)

        !----- alltoall communication of sendbuf to recvbuf
        ierr = cudaStreamSynchronize()
        call MPI_Alltoall(p%sendbuf, nblksize, MPI_DOUBLE_PRECISION, &
                          p%recvbuf, nblksize, MPI_DOUBLE_PRECISION, &
                          p%ptdma_world, ierr)

        call mem_unite_slab_xy<<<p%blocks_alltoall, p%threads>>> &
                                (p%recvbuf, slab_xy, p%nx_sys, p%ny_sys, p%nz_row_rd, p%nprocs)

    end subroutine transpose_slab_yz_to_xy

    !>
    !> @brief   Subroutine to rearrange y-z slab to a 1D array for MPI_Alltoall
    !> @param   slab_yz     Coefficient array in the shape of y-z slab
    !> @param   array1D     1-D Coefficient array
    !> @param   n1, n2, n3  Dimension of array 'slab_xy'
    !> @param   nprocs      Number of processes in 'mpi_comm'
    !>
    attributes(global) subroutine mem_detach_slab_yz(slab_yz, array1D, n1, n2, n3, nprocs)

        implicit none

        integer, value, intent(in)  :: n1, n2, n3, nprocs
        double precision, device, intent(in)    :: slab_yz(:, :, :)
        double precision, device, intent(out)   :: array1D(:)

        ! Variables to calculate indices for in and out arrays
        integer :: i, j, k, kblk, blksize
        integer :: pos_k, pos_i, pos_j, pos

        i = (blockidx%x - 1) * blockdim%x + threadidx%x
        j = (blockidx%y - 1) * blockdim%y + threadidx%y
        k = (blockidx%z - 1) * blockdim%z + threadidx%z

        pos_i = (i - 1) * n2 * n3
        pos_j = (j - 1) * n3
        blksize  = n1 * n2 * n3 / nprocs

        do kblk = 1, nprocs
            pos_k = k + (kblk - 1) * blksize
            pos   = pos_k + pos_j + pos_i
            array1D(pos) = slab_yz(i, j, k + (kblk - 1) * n3)
        enddo

    end subroutine mem_detach_slab_yz

    !>
    !> @brief   Subroutine to rearrange a 1D array to x-y slab after MPI_Alltoall
    !> @param   array1D     1-D Coefficient array
    !> @param   slab_xy     Coefficient array in the shape of x-y slab
    !> @param   n1, n2, n3  Dimension of array 'slab_xy'
    !> @param   nprocs      Number of processes in 'mpi_comm'
    !>
    attributes(global) subroutine mem_unite_slab_xy(array1D, slab_xy, n1, n2, n3, nprocs)

        implicit none

        integer, value, intent(in)  :: n1, n2, n3, nprocs
        double precision, device, intent(in)    :: array1D(:)
        double precision, device, intent(out)   :: slab_xy(:, :, :)

        ! Variables to calculate indices for in and out arrays
        integer :: i, j, k, kblk, n1blksize, blksize
        integer :: pos, pos_i, pos_j, pos_k

        n1blksize  = n1 / nprocs
        blksize  = n1blksize * n2 * n3

        i = (blockidx%x - 1) * blockdim%x + threadidx%x
        j = (blockidx%y - 1) * blockdim%y + threadidx%y
        k = (blockidx%z - 1) * blockdim%z + threadidx%z

        pos_i = (i - 1) * n2 * n3
        pos_j = (j - 1) * n3

        do kblk = 1, nprocs
            pos_k = k + (kblk - 1) * blksize
            pos   = pos_k + pos_j + pos_i
            slab_xy(i + (kblk - 1) * n1blksize, j, k) = array1D(pos)
        enddo

    end subroutine mem_unite_slab_xy

end module PaScaL_TDMA_cuda