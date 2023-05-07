!======================================================================================================================
!> @file        tdmas.f90
!> @brief       Tridiagonal matrix (TDM) solvers using the Thomas algorithm.
!> @details     A single TDM solver and many TDM solver with non-cyclic and cyclic conditions.
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
!> @brief       Solve a single tridiagonal system of equations using the Thomas algorithm.
!> @param       a       Coefficients in lower diagonal elements
!> @param       b       Coefficients in diagonal elements
!> @param       c       Coefficients in upper diagonal elements
!> @param       d       Coefficients in the right-hand side terms
!> @param       n1      Number of rows in each process, dimension of tridiagonal matrix N divided by nprocs
!>
subroutine tdma_single(a, b, c, d, n1)

    implicit none

    integer, intent(in) :: n1
    double precision, intent(inout) :: a(n1), b(n1), c(n1), d(n1)
    
    integer :: i
    double precision :: r

    d(1)=d(1)/b(1)
    c(1)=c(1)/b(1)

    do i=2,n1
        r=1.d0/(b(i)-a(i)*c(i-1))
        d(i)=r*(d(i)-a(i)*d(i-1))
        c(i)=r*c(i)
    enddo

    do i=n1-1,1,-1
        d(i)=d(i)-c(i)*d(i+1)
    enddo

end subroutine tdma_single

!>
!> @brief       Solve a single cyclic tridiagonal system of equations using the Thomas algorithm.
!> @param       a       Coefficients in lower diagonal elements
!> @param       b       Coefficients in diagonal elements
!> @param       c       Coefficients in upper diagonal elements
!> @param       d       Coefficients in the right-hand side terms
!> @param       n1      Number of rows in each process, dimension of tridiagonal matrix N divided by nprocs
!>
subroutine tdma_cycl_single(a, b, c, d, n1)

    implicit none

    integer, intent(in) :: n1
    double precision, intent(inout) :: a(n1), b(n1), c(n1), d(n1)
    
    integer :: i
    double precision :: rr
    double precision, allocatable, dimension(:) :: e

    allocate(e(1:n1))

    e(:)=0
    e(2 ) = -a(2 )
    e(n1) = -c(n1)

    d(2)=d(2)/b(2)
    e(2)=e(2)/b(2)
    c(2)=c(2)/b(2)

    do i=3,n1
        rr=1.d0/(b(i)-a(i)*c(i-1))
        d(i)=rr*(d(i)-a(i)*d(i-1))
        e(i)=rr*(e(i)-a(i)*e(i-1))
        c(i)=rr*c(i)
    enddo

    do i=n1-1,2,-1
        d(i)=d(i)-c(i)*d(i+1)
        e(i)=e(i)-c(i)*e(i+1)
    enddo

    d(1)=(d(1)-a(1)*d(n1)-c(1)*d(2))/(b(1)+a(1)*e(n1)+c(1)*e(2))

    do i = 2, n1
        d(i) = d(i) + d(1)*e(i)
    end do

    deallocate(e)
end subroutine tdma_cycl_single

!>
!> @brief       Solve many tridiagonal systems of equations using the Thomas algorithm.
!>              First index indicates the number of independent many tridiagonal systems to use vectorization.
!>              Second index indicates the row number in the tridiagonal system .
!> @param       a       Coefficient array in lower diagonal elements
!> @param       b       Coefficient array in diagonal elements
!> @param       c       Coefficient array in upper diagonal elements
!> @param       d       Coefficient array in the right-hand side terms
!> @param       n1      Number of tridiagonal systems per process
!> @param       n2      Number of rows in each process, size of the tridiagonal matrix N divided by nprocs
!>
subroutine tdma_many(a, b, c, d, n1, n2)

    implicit none

    integer, intent(in) :: n1,n2
    double precision, intent(inout) :: a(n1,n2), b(n1,n2), c(n1,n2), d(n1,n2)
    
    integer :: i,j
    double precision, allocatable, dimension(:) :: r

    allocate(r(1:n1))

    do i=1,n1
        d(i,1)=d(i,1)/b(i,1)
        c(i,1)=c(i,1)/b(i,1)
    enddo

    do j=2,n2
        do i=1,n1
            r(i)=1.d0/(b(i,j)-a(i,j)*c(i,j-1))
            d(i,j)=r(i)*(d(i,j)-a(i,j)*d(i,j-1))
            c(i,j)=r(i)*c(i,j)
        enddo
    enddo

    do j=n2-1,1,-1
        do i=1,n1
            d(i,j)=d(i,j)-c(i,j)*d(i,j+1)
        enddo
    enddo

    deallocate(r)

end subroutine tdma_many

!>
!> @brief       Solve many cyclic tridiagonal systems of equations using the Thomas algorithm.
!>              First index indicates the number of independent many tridiagonal systems to use vectorization.
!>              Second index indicates the row number in the tridiagonal system.
!> @param       a       Coefficient array in lower diagonal elements
!> @param       b       Coefficient array in diagonal elements
!> @param       c       Coefficient array in upper diagonal elements
!> @param       d       Coefficient array in the right-hand side terms
!> @param       n1      Number of tridiagonal systems per process
!> @param       n2      Number of rows in each process, size of the tridiagonal matrix N divided by nprocs
!>
subroutine tdma_cycl_many(a, b, c, d, n1, n2)

    implicit none

    integer, intent(in) :: n1,n2
    double precision, intent(inout) :: a(n1,n2), b(n1,n2), c(n1,n2), d(n1,n2)
    
    integer :: i,j
    double precision, allocatable, dimension(:,:) :: e
    double precision, allocatable, dimension(:) :: rr

    allocate(e(1:n1,1:n2),rr(1:n1))

    do i=1,n1
        e(i,:)=0
        e(i,2 ) = -a(i,2 )
        e(i,n2) = -c(i,n2)
    enddo

    do i=1,n1
        d(i,2)=d(i,2)/b(i,2)
        e(i,2)=e(i,2)/b(i,2)
        c(i,2)=c(i,2)/b(i,2)
    enddo

    do j=3,n2
        do i=1,n1
            rr(i)=1.d0/(b(i,j)-a(i,j)*c(i,j-1))
            d(i,j)=rr(i)*(d(i,j)-a(i,j)*d(i,j-1))
            e(i,j)=rr(i)*(e(i,j)-a(i,j)*e(i,j-1))
            c(i,j)=rr(i)*c(i,j)
        enddo
    enddo

    do j=n2-1,2,-1
        do i=1,n1
            d(i,j)=d(i,j)-c(i,j)*d(i,j+1)
            e(i,j)=e(i,j)-c(i,j)*e(i,j+1)
        enddo
    enddo

    do i=1,n1
        d(i,1)=(d(i,1)-a(i,1)*d(i,n2)-c(i,1)*d(i,2))/(b(i,1)+a(i,1)*e(i,n2)+c(i,1)*e(i,2))
    enddo

    do j = 2,n2
        do i=1,n1
            d(i,j) = d(i,j) + d(i,1)*e(i,j)
        enddo
    end do
    
    deallocate(e,rr)

end subroutine tdma_cycl_many
