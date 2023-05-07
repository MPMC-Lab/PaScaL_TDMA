!=======================================================================================================================
!> @file        tdmas_cuda.cuf
!> @brief       CUDA version of tridiagonal matrix (TDM) solvers using the Thomas algorithm.
!> @details     A single TDM solver and many TDM solver with non-cyclic and cyclic conditions.
!>
!> @author      
!>              - Mingyu Yang (yang926@yonsei.ac.kr), Multi-Physics Modeling and Computation Lab., Yonsei University
!>              - Ji-Hoon Kang (jhkang@kisti.re.kr), Korea Institute of Science and Technology Information
!>              - Ki-Ha Kim (k-kiha@yonsei.ac.kr), Department of Computational Science & Engineering, Yonsei University
!>              - Jung-Il Choi (jic@yonsei.ac.kr), Department of Computational Science & Engineering, Yonsei University
!>
!> @date        March 2023
!> @version     2.0
!> @par         Copyright
!>              Copyright (c) 2019-2023 Mingyu Yang, Ki-Ha Kim and Jung-Il choi, Yonsei University and 
!>              Ji-Hoon Kang, Korea Institute of Science and Technology Information, All rights reserved.
!> @par         License     
!>              This project is release under the terms of the MIT License (see LICENSE in )
!=======================================================================================================================

module tdma

contains

    !>
    !> @brief   Solve many tridiagonal systems of equations using the Thomas algorithm.
    !>          First & second indices indicate the number of independent many tridiagonal systems for parallelization.
    !>          Third index indicates the row number in the tridiagonal system.
    !> @param   a           Coefficient array in lower diagonal elements
    !> @param   b           Coefficient array in diagonal elements
    !> @param   c           Coefficient array in upper diagonal elements
    !> @param   d           Coefficient array in the right-hand side terms
    !> @param   nz_row      Row size of tridiagonal matrix in z-direction that is the solving direction
    !>
    attributes(global) subroutine tdma_many_cuda(a, b, c, d, nz_row)

        implicit none

        integer, value, intent(in)      :: nz_row
        double precision, device, intent(in)    :: a(:, :, :), b(:, :, :)
        double precision, device, intent(inout) :: c(:, :, :), d(:, :, :)
        
        integer :: i, j, k
        integer :: tj, tk
        double precision :: r

        double precision, shared :: a1(blockdim%x + 1, blockdim%y)
        double precision, shared :: b1(blockdim%x + 1, blockdim%y)
        double precision, shared :: c0(blockdim%x + 1, blockdim%y), c1(blockdim%x + 1, blockdim%y)
        double precision, shared :: d0(blockdim%x + 1, blockdim%y), d1(blockdim%x + 1, blockdim%y)

        ! Global index
        j = (blockidx%x - 1) * blockdim%x + threadidx%x
        k = (blockidx%y - 1) * blockdim%y + threadidx%y

        ! Local index in block
        tj = threadidx%x
        tk = threadidx%y

        ! Using shared memory
		! First and second indices are for thread IDs
        b1(tj, tk) = b(j, k, 1)
        c1(tj, tk) = c(j, k, 1)
        d1(tj, tk) = d(j, k, 1)

        d1(tj, tk) = d1(tj, tk) / b1(tj, tk)
        c1(tj, tk) = c1(tj, tk) / b1(tj, tk)
 
        d(j,k,1) = d1(tj, tk)
        c(j,k,1) = c1(tj, tk)

        do i = 2, nz_row

            ! Pipeline copy of (i-1)th data using shared memory
            c0(tj, tk) = c1(tj, tk)
            d0(tj, tk) = d1(tj, tk)

            ! Load i-th data from global memory
            a1(tj, tk) = a(j, k, i)
            b1(tj, tk) = b(j, k, i)
            c1(tj, tk) = c(j, k, i)
            d1(tj, tk) = d(j, k, i)

            r = 1.0d0 / ( b1(tj, tk) - a1(tj, tk) * c0(tj, tk))
            d1(tj, tk) = r * (d1(tj, tk) - a1(tj, tk) * d0(tj, tk))
            c1(tj, tk) = r * c1(tj, tk)

            ! Save updated i-th data to global memory
            d(j,k,i) = d1(tj,tk)
            c(j,k,i) = c1(tj,tk)

        enddo

        do i = nz_row - 1, 1, -1
            c0(tj, tk) = c(j, k, i)
            d0(tj, tk) = d(j, k, i)
            d0(tj, tk) = d0(tj, tk) - c0(tj, tk) * d1(tj, tk)
            d1(tj, tk) = d0(tj, tk)
            d(j, k, i) = d0(tj, tk)
        enddo

    end subroutine tdma_many_cuda

    !>
    !> @brief   Solve many tridiagonal systems of equations using the Thomas algorithm.
    !>          First & second indices indicate the number of independent many tridiagonal systems for parallelization.
    !>          Third index indicates the row number in the tridiagonal system.
    !> @param   a           Coefficient array in lower diagonal elements
    !> @param   b           Coefficient array in diagonal elements
    !> @param   c           Coefficient array in upper diagonal elements
    !> @param   d           Coefficient array in the right-hand side terms
    !> @param   nz_row      Row size of tridiagonal matrix in z-direction that is the solving direction
    !>
    attributes(global) subroutine tdma_cycl_many_cuda(a, b, c, d, e, nz_row)

        implicit none

        integer, value, intent(in)      :: nz_row
        double precision, device, intent(in)    :: a(:, :, :), b(:, :, :)
        double precision, device, intent(inout) :: c(:, :, :), d(:, :, :), e(:, :, :)

        integer :: i, j, k
        integer :: tj, tk
        double precision :: r

        double precision, shared :: a1(blockdim%x + 1, blockdim%y)
        double precision, shared :: b1(blockdim%x + 1, blockdim%y)
        double precision, shared :: c0(blockdim%x + 1, blockdim%y), c1(blockdim%x + 1, blockdim%y)
        double precision, shared :: d0(blockdim%x + 1, blockdim%y), d1(blockdim%x + 1, blockdim%y)
        double precision, shared :: e0(blockdim%x + 1, blockdim%y), e1(blockdim%x + 1, blockdim%y)

        ! Global index
        j = (blockidx%x - 1) * blockdim%x + threadidx%x
        k = (blockidx%y - 1) * blockdim%y + threadidx%y

        ! Local index in block
        tj = threadidx%x
        tk = threadidx%y

        do i = 1, nz_row
            e(j, k, i)  = 0.0d0
        enddo

        ! Using shared memory
		! First and second indices are for thread IDs
        e(j, k, 2)  = -a(j, k, 2)
        e(j, k, nz_row) = -c(j, k, nz_row)

        d1(tj, tk) = d(j, k, 2)
        b1(tj, tk) = b(j, k, 2)
        c1(tj, tk) = c(j, k, 2)
        e1(tj, tk) = e(j, k, 2)

        d1(tj, tk) = d1(tj, tk) / b1(tj, tk)
        c1(tj, tk) = c1(tj, tk) / b1(tj, tk)
        e1(tj, tk) = e1(tj, tk) / b1(tj, tk)

        d(j, k, 2)  = d1(tj, tk)
        c(j, k, 2)  = c1(tj, tk)
        e(j, k, 2)  = e1(tj, tk)

        do i = 3, nz_row
            ! Pipeline copy of (i-1)th data using shared memory
            c0(tj, tk) = c1(tj, tk)
            d0(tj, tk) = d1(tj, tk)
            e0(tj, tk) = e1(tj, tk)

            ! Load i-th data from global memory
            a1(tj, tk) = a(j, k, i)
            b1(tj, tk) = b(j, k, i)
            c1(tj, tk) = c(j, k, i)
            d1(tj, tk) = d(j, k, i)
            e1(tj, tk) = e(j, k, i)

            r = 1.0d0 / (b1(tj, tk) - a1(tj, tk) * c0(tj, tk))
            d1(tj, tk)= r * (d1(tj, tk) - a1(tj, tk) * d0(tj, tk))
            e1(tj, tk)= r * (e1(tj, tk) - a1(tj, tk) * e0(tj, tk))
            c1(tj, tk)= r * c1(tj, tk)

            ! Save updated i-th data to global memory
            e(j,k, i) = e1(tj, tk)
            d(j,k, i) = d1(tj, tk)
            c(j,k, i) = c1(tj, tk)
        enddo

        do i = nz_row-1, 2, -1
            ! Load i-th data from global memory
            c0(tj, tk) = c(j, k, i)
            d0(tj, tk) = d(j, k, i)
            e0(tj, tk) = e(j, k, i)

            d0(tj, tk) = d0(tj, tk) - c0(tj, tk) * d1(tj, tk)
            e0(tj, tk) = e0(tj, tk) - c0(tj, tk) * e1(tj, tk) 

            ! Pipeline copy of updated i-th data using shared memory
            d1(tj, tk) = d0(tj, tk)
            e1(tj, tk) = e0(tj, tk)
            d(j, k, i) = d0(tj, tk)
            e(j, k, i) = e0(tj, tk)
        enddo

        a1(tj, tk) = a(j, k, 1)
        b1(tj, tk) = b(j, k, 1)
        c1(tj, tk) = c(j, k, 1)
        d1(tj, tk) = d(j, k, 1)
        e1(tj, tk) = e(j, k, nz_row)

        d1(tj, tk) = (d1(tj, tk) - a1(tj, tk) * d(j, k, nz_row) - c1(tj, tk) * d0(tj, tk)) / &
                     (b1(tj, tk) + a1(tj, tk) * e1(tj, tk) + c1(tj, tk) * e0(tj, tk))
        d(j,k, 1) = d1(tj, tk)

        do i = 2, nz_row
            e1(tj, tk) = e(j, k, i)
            d(j,k, i) = d(j, k, i) + d1(tj, tk) * e1(tj, tk)
        enddo

    end subroutine tdma_cycl_many_cuda

end module tdma