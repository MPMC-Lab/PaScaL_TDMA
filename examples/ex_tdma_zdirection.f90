program main
    use mpi                        ! MPI parallel environment
    use cudafor                    ! CUDA Fortran interfaces
    use PaScaL_TDMA_cuda           ! PaScaL-TDMA CUDA solver module
    implicit none

    ! MPI & CUDA process information
    integer :: ierr, myrank, nprocs
    integer :: gpurank, ngpu

    ! Global domain sizes
    integer :: n1, n2, n3
    ! Local subdomain sizes (assigned by domain decomposition)
    integer :: n1sub, n2sub, n3sub

    ! Host arrays for tridiagonal system (Aa, Ab, Ac) and RHS (B)
    real*8, allocatable, dimension(:,:,:) :: Aa,Ab,Ac,B
    ! Device (GPU) arrays for the same data
    real*8, allocatable, dimension(:,:,:), device :: Aa_d, Ab_d, Ac_d, B_d

    ! PaScaL-TDMA plan structure (contains buffers, streams, communicator info)
    type(ptdma_plan_cuda) :: exampleplan
    
    ! Thread settings for modified-Thomas stage and reduced system stage
    integer :: nthread_modithomas, nthread_reduced

    ! z-direction decomposition indices for each MPI rank
    integer :: ia,ib

    !===========================================================
    ! MPI INIT
    !===========================================================
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    if(myrank==0) write(*,*) " MPI processes:", nprocs
    ! Initializes MPI and retrieves the number of processes and the local rank.
    ! All MPI processes will each solve a slab (z-range) of the global domain.

    !===========================================================
    ! CUDA INIT
    !===========================================================
    ierr = CUDAGETDEVICECOUNT(ngpu)      ! Query number of available GPUs
    gpurank = mod(myrank,ngpu)           ! Map MPI rank → GPU index
    ierr = CUDASETDEVICE(gpurank)        ! Bind MPI rank to its GPU
    ierr = CUDADEVICESYNCHRONIZE()       ! Ensure device is ready
    if(myrank==0) write(*,*) " CUDA devices :", ngpu
    ! Each MPI rank selects a GPU device. Important for multi-GPU clusters.

    !===========================================================
    ! Set grid sizes
    !===========================================================
    n1 = 64                               ! x-dimension
    n2 = 64                               ! y-dimension
    n3 = 2048                             ! z-dimension (split across MPI ranks)

    ! CUDA thread configuration for PaScaL-TDMA
    nthread_modithomas = 128
    nthread_reduced    = 128
    if(myrank==0) write(*,*) " Grid size    :", n1,n2,n3

    !===========================================================
    ! Domain decomposition
    !===========================================================
    call para(0,n3-1,nprocs,myrank,ia,ib) ! Compute local z-range [ia, ib]
    n1sub = n1
    n2sub = n2
    n3sub = ib - ia + 1                   ! Local domain thickness
    if(myrank==0) write(*,*) " Subdomain z-size:", n3sub
    ! Each rank gets a contiguous z-slice of the full 3D domain.

    !===========================================================
    ! Allocate memory
    !===========================================================
    allocate(Aa(0:n1sub-1,0:n2sub-1,0:n3sub-1), Aa_d(0:n1sub-1,0:n2sub-1,0:n3sub-1))
    allocate(Ab(0:n1sub-1,0:n2sub-1,0:n3sub-1), Ab_d(0:n1sub-1,0:n2sub-1,0:n3sub-1))
    allocate(Ac(0:n1sub-1,0:n2sub-1,0:n3sub-1), Ac_d(0:n1sub-1,0:n2sub-1,0:n3sub-1))
    allocate(B (0:n1sub-1,0:n2sub-1,0:n3sub-1), B_d (0:n1sub-1,0:n2sub-1,0:n3sub-1))
    if(myrank==0) write(*,*) " Memory allocated"
    ! Host and device memory for the tridiagonal system.
    ! Arrays follow CSR-like structure but stored as 3D blocks.

    !===========================================================
    ! Initialize matrix system
    !===========================================================
    Aa(:,:,:) = dble( 1)       ! Lower diagonal coefficient
    Ab(:,:,:) = dble(-2)       ! Main diagonal
    Ac(:,:,:) = dble( 1)       ! Upper diagonal
    B (:,:,:) = dble( 0)       ! RHS vector

    ! Boundary conditions at the physical domain edges
    if(myrank==0       ) B(:,:,0      ) = dble(-1)
    if(myrank==nprocs-1) B(:,:,n3sub-1) = dble(-1)
    if(myrank==0) write(*,*) " Matrix system initialized"

    ! Copy host data → device
    Aa_d = Aa
    Ab_d = Ab
    Ac_d = Ac
    B_d  = B
    if(myrank==0) write(*,*) " Data copied to device"
    ! At this point, GPU memory contains the full local tridiagonal system.

    !===========================================================
    ! *** Pascal TDMA example ***
    !===========================================================
    if(myrank==0) write(*,*) " Starting PascaL TDMA solver"

    ! Create solver plan (communication buffers, streams, workspaces)
    call pascal_plan_create(exampleplan, (n1sub*n2sub), MPI_COMM_WORLD, myrank, nprocs, &
                            nthread_modithomas, nthread_reduced)

    ! Execute multi-GPU distributed tridiagonal solve along z-direction
    call pascal_solver(exampleplan, Aa_d, Ab_d, Ac_d, B_d, (n1sub*n2sub), n3sub)

    ! Release GPU buffers and internal structures
    call pascal_plan_clean(exampleplan)

    if(myrank==0) write(*,*) " PascaL TDMA solver finished"

    !===========================================================
    ! Print small portion of solution
    !===========================================================
    B = B_d              ! Copy device → host for printing
    call checkprint(B, n1sub, n2sub, n3sub, ia, myrank, nprocs)
    ! checkprint() prints a few lines from the solution for debugging.

    !===========================================================
    ! FINALIZE
    !===========================================================
    deallocate(Aa,Aa_d)
    deallocate(Ab,Ab_d)
    deallocate(Ac,Ac_d)
    deallocate(B ,B_d)
    call MPI_FINALIZE(ierr)
    ! Clean exit for both MPI and CUDA memory.
end


!=======================================================================
! checkprint : prints selected entries from the solution array B.
! This function is primarily for verification and debugging.
! Each MPI process prints its own assigned slab in rank order.
!=======================================================================
subroutine checkprint(B, n1sub, n2sub, n3sub, ia, myrank, nprocs)
    use mpi
    implicit none

    ! Arguments
    real*8, dimension(0:n1sub-1,0:n2sub-1,0:n3sub-1) :: B
    integer :: n1sub, n2sub, n3sub
    integer :: ia
    integer :: myrank, nprocs

    ! Local variables
    integer :: iter, k, ierr

    do iter = 0, nprocs-1
        if(myrank == iter) then
            write(*,'(1A5,1I5,1A2)') " Rank:", myrank, "--"

            ! Print first few z-indices
            do k = 0, 2
                write(*,'(1I8,1I8,F8.3,A4,F8.3,A4,F8.3)')                     &
                     k, k+ia, B(0,0,k), " ...", B(n1sub/2,n2sub/2,k), " ...", &
                     B(n1sub-1,n2sub-1,k)
            end do

            write(*,'(1A8)') "     ..."

            ! Print last two z-indices
            do k = n3sub-2, n3sub-1
                write(*,'(1I8,1I8,F8.3,A4,F8.3,A4,F8.3)')                     &
                     k, k+ia, B(0,0,k), " ...", B(n1sub/2,n2sub/2,k), " ...", &
                     B(n1sub-1,n2sub-1,k)
            end do
        endif

        ! Synchronize between MPI processes for ordered printing
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    end do

end subroutine checkprint
