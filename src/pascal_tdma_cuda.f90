module PaScaL_TDMA_cuda
    use mpi                      ! MPI for distributed communication
    use cudafor                  ! CUDA Fortran extensions (device, dim3, etc.)
    implicit none

    !===================================================================
    ! ptdma_plan_cuda
    ! ----------------
    ! This derived type stores all metadata and device buffers required
    ! by the parallel TDMA solver:
    !   - MPI communicator and rank layout
    !   - Problem sizes for reduced (rd) and transformed (tr) systems
    !   - All gather / scatter descriptors for communication
    !   - Packed communication buffers (host/device)
    !   - CUDA launch configurations (block / grid)
    !   - Device work arrays for reduced and transformed systems
    !===================================================================
    type, public :: ptdma_plan_cuda
        ! MPI context
        integer :: ptdma_world,myrank,nprocs

        ! Global and local sizes of reduced system (Ard,Brd,Crd,Drd)
        ! Nrd_global = (/ Nsys, 2 /)
        ! Nrd_local  = (/ local_Nsys, 2 /)
        integer :: Nrd_global(0:1),Nrd_local(0:1)

        ! Global and local sizes of transformed (reduced) system (Atr,Btr,Ctr,Dtr)
        ! Ntr_global = (/ local_Nsys, 2*nprocs /)
        ! Ntr_local  = (/ local_Nsys, 2 /)
        integer :: Ntr_global(0:1),Ntr_local(0:1)

        ! Host-side gather descriptors:
        !   gather_N*_local(0,:)  : local "row" (first dimension) size
        !   gather_N*_local(1,:)  : local "column" (second dimension) size
        !   gather_N*_start(0,:)  : starting row index for each rank
        !   gather_N*_start(1,:)  : starting column index for each rank
        integer, allocatable, dimension(:,:) :: gather_Nrd_local,gather_Ntr_local
        integer, allocatable, dimension(:,:) :: gather_Nrd_start,gather_Ntr_start

        ! Device-side copies of the same gather descriptors
        integer, allocatable, dimension(:,:), device :: gather_Nrd_local_d,gather_Ntr_local_d
        integer, allocatable, dimension(:,:), device :: gather_Nrd_start_d,gather_Ntr_start_d

        ! Host-side buffer size and displacement descriptors for A-side
        integer, allocatable, dimension(:) :: bufsubsize_A,bufstart_A
        integer, allocatable, dimension(:) :: BIGbufsubsize_A,BIGbufstart_A

        ! Host-side buffer size and displacement descriptors for B-side
        integer, allocatable, dimension(:) :: bufsubsize_B,bufstart_B
        integer, allocatable, dimension(:) :: BIGbufsubsize_B,BIGbufstart_B

        ! Device-side large communication buffers for packing/unpacking
        real*8, allocatable, dimension(:), device :: BIGbuf_A, BIGbuf_B 

        ! CUDA launch configuration for:
        !   t_* : thread block dimensions (dim3)
        !   b_* : grid dimensions       (dim3)
        !   tdma     : local TDMA sweeps
        !   rdtdma   : reduced TDMA sweeps
        !   pack     : pack data into communication buffers
        !   unpack   : unpack data from communication buffers
        !   intiAb   : initialize reduced system (e.g. identity / ones)
        type(dim3) :: t_tdma, t_rdtdma, t_pack, t_unpack, t_intiAb
        type(dim3) :: b_tdma, b_rdtdma, b_pack, b_unpack, b_intiAb

        ! Device-side reduced system (size Nsys x 2):
        !   Ard, Brd, Crd, Drd represent boundary-compressed TDMA coefficients
        real*8, allocatable, dimension(:,:), device :: Ard,Brd,Crd,Drd

        ! Device-side transformed system (size local_N x 2*nprocs):
        !   Atr, Btr, Ctr, Dtr form the globally reduced tridiagonal system
        real*8, allocatable, dimension(:,:), device :: Atr,Btr,Ctr,Dtr
        
    end type ptdma_plan_cuda

    !===================================================================
    ! ptdma_timing_cuda
    ! -----------------
    ! Phase timing container for the profiled solver entry point.
    ! Times are local to each MPI rank. Study drivers reduce them across
    ! ranks when writing CSV output.
    !===================================================================
    type, public :: ptdma_timing_cuda
        real(8) :: total = 0.0d0
        real(8) :: local_compute = 0.0d0
        real(8) :: pack_forward = 0.0d0
        real(8) :: mpi_forward = 0.0d0
        real(8) :: unpack_forward = 0.0d0
        real(8) :: reduced_compute = 0.0d0
        real(8) :: pack_backward = 0.0d0
        real(8) :: mpi_backward = 0.0d0
        real(8) :: unpack_backward = 0.0d0
        real(8) :: update_compute = 0.0d0
    end type ptdma_timing_cuda
    
contains

    subroutine pascal_timing_reset(timing)
        implicit none
        type(ptdma_timing_cuda) :: timing

        timing%total = 0.0d0
        timing%local_compute = 0.0d0
        timing%pack_forward = 0.0d0
        timing%mpi_forward = 0.0d0
        timing%unpack_forward = 0.0d0
        timing%reduced_compute = 0.0d0
        timing%pack_backward = 0.0d0
        timing%mpi_backward = 0.0d0
        timing%unpack_backward = 0.0d0
        timing%update_compute = 0.0d0
    end subroutine pascal_timing_reset

    !===================================================================
    ! pascal_plan_create
    ! -------------------
    ! Set up the ptdma_plan_cuda structure:
    !   - Compute local domain size using "para"
    !   - Allocate device work arrays (Ard,Brd,Crd,Drd,Atr,Btr,Ctr,Dtr)
    !   - Build gather descriptors for reduced and transformed systems
    !   - Construct host/device communication buffers and offsets
    !   - Determine CUDA kernel launch configurations (t_*, b_*)
    !
    ! Input:
    !   plan       : solver plan to be initialized
    !   Nsys       : number of independent tridiagonal systems
    !   commworld  : MPI communicator
    !   myrank     : MPI rank of this process
    !   nprocs     : total number of MPI processes
    !   tmp_opti1  : preferred thread block size for local TDMA
    !   tmp_opti2  : preferred thread block size for reduced TDMA
    !===================================================================
    subroutine pascal_plan_create(plan,Nsys,commworld,myrank,nprocs,tmp_opti1,tmp_opti2)
        implicit none
        type(ptdma_plan_cuda) :: plan 
        integer :: Nsys
        integer :: commworld,myrank,nprocs
        integer :: i,ia,ib

        integer :: tmp_N,tmp_opti1,tmp_opti2,tmp_opti3
        integer, allocatable, dimension(:) :: tmp_int
        allocate(tmp_int(0:nprocs-1))

        ! Store MPI context inside the plan
        plan%ptdma_world = commworld
        plan%myrank = myrank
        plan%nprocs = nprocs

        ! Compute local range of systems [ia, ib] for this rank
        call para(0,Nsys-1,nprocs,myrank,ia,ib)
        tmp_N = ib-ia+1

        ! Global and local grid sizes for reduced and transformed systems
        plan%Nrd_global(0:1) = (/Nsys,2/)
        plan%Ntr_global(0:1) = (/tmp_N,2*nprocs/)

        plan%Nrd_local(0:1)  = (/tmp_N,2/)
        plan%Ntr_local(0:1)  = (/tmp_N,2/)

        ! Device arrays:
        !   Ard,Brd,Crd,Drd : reduced coefficients (Nsys x 2)
        allocate(plan%Ard(0:Nsys-1,0:1),plan%Brd(0:Nsys-1,0:1),plan%Crd(0:Nsys-1,0:1),plan%Drd(0:Nsys-1,0:1))

        !   Atr,Btr,Ctr,Dtr : transformed system over Nprocs (tmp_N x 2*nprocs)
        allocate(plan%Atr(0:tmp_N-1,0:2*nprocs-1),plan%Btr(0:tmp_N-1,0:2*nprocs-1))
        allocate(plan%Ctr(0:tmp_N-1,0:2*nprocs-1),plan%Dtr(0:tmp_N-1,0:2*nprocs-1))
        
        ! Host-side gather information (local sizes per rank)
        allocate(plan%gather_Nrd_local(0:1,0:nprocs-1),plan%gather_Ntr_local(0:1,0:nprocs-1))
        allocate(plan%gather_Nrd_local_d(0:1,0:nprocs-1),plan%gather_Ntr_local_d(0:1,0:nprocs-1))

        ! Host-side gather start indices (prefix sums for packed layout)
        allocate(plan%gather_Nrd_start(0:1,0:nprocs-1),plan%gather_Ntr_start(0:1,0:nprocs-1))
        allocate(plan%gather_Nrd_start_d(0:1,0:nprocs-1),plan%gather_Ntr_start_d(0:1,0:nprocs-1))

        !----------------------------------------------------------------
        ! Build per-rank descriptors for reduced and transformed systems
        !----------------------------------------------------------------
        do i = 0, nprocs-1
            ! Local problem size for rank i
            call para(0,Nsys-1,nprocs,i,ia,ib)
            plan%gather_Nrd_local(0:1,i) = (/ib-ia+1,2/)
            plan%gather_Ntr_local(0:1,i) = plan%Ntr_local(0:1)

            ! Compute starting indices (prefix sums) for each rank
            plan%gather_Nrd_start(0:1,i) = (/sum(plan%gather_Nrd_local(0,0:i))-plan%gather_Nrd_local(0,i),0/) 
            plan%gather_Ntr_start(0:1,i) = (/0,sum(plan%gather_Ntr_local(1,0:i))-plan%gather_Ntr_local(1,i)/)

            tmp_int(i) = ib-ia+1
        end do
        
        ! Copy gather descriptors to device
        plan%gather_Nrd_local_d= plan%gather_Nrd_local
        plan%gather_Ntr_local_d= plan%gather_Ntr_local
        plan%gather_Nrd_start_d= plan%gather_Nrd_start
        plan%gather_Ntr_start_d= plan%gather_Ntr_start
        
        !----------------------------------------------------------------
        ! Allocate host-side buffer descriptors for all-to-all communication
        !----------------------------------------------------------------
        allocate(plan%bufsubsize_A(0:nprocs-1),plan%bufsubsize_B(0:nprocs-1))
        allocate(plan%bufstart_A(0:nprocs-1),plan%bufstart_B(0:nprocs-1))
        allocate(plan%BIGbufsubsize_A(0:nprocs-1),plan%BIGbufsubsize_B(0:nprocs-1))
        allocate(plan%BIGbufstart_A(0:nprocs-1),plan%BIGbufstart_B(0:nprocs-1))

        ! Per-rank buffer extents and starts for reduced (A) and transformed (B) data
        do i = 0, nprocs-1
            plan%bufsubsize_A(i) = plan%gather_Nrd_local(0,i)*plan%gather_Nrd_local(1,i)
            plan%bufsubsize_B(i) = plan%gather_Ntr_local(0,i)*plan%gather_Ntr_local(1,i)

            plan%bufstart_A(i) = sum(plan%bufsubsize_A(0:i)) - plan%bufsubsize_A(i)
            plan%bufstart_B(i) = sum(plan%bufsubsize_B(0:i)) - plan%bufsubsize_B(i)
        end do

        ! "BIG" buffers are sized at 3× the minimal requirement to store
        ! three coefficient arrays (e.g. A, C, D) per rank.
        plan%BIGbufsubsize_A(:) = 3*plan%bufsubsize_A
        plan%BIGbufsubsize_B(:) = 3*plan%bufsubsize_B
        plan%BIGbufstart_A(:)   = 3*plan%bufstart_A
        plan%BIGbufstart_B(:)   = 3*plan%bufstart_B

        ! Device-side consolidated communication buffers
        allocate(plan%BIGbuf_A(0:sum(plan%BIGbufsubsize_A(:))-1),plan%BIGbuf_B(0:sum(plan%BIGbufsubsize_B(:))-1))


        !----------------------------------------------------------------
        ! Configure CUDA thread block sizes
        !----------------------------------------------------------------
        ! tmp_opti1 = Nsys
        ! tmp_opti2 = tmp_N
        tmp_opti3 = maxval(tmp_int(:))
        ! if(Nsys>512) tmp_opti1 = 512
        ! if(tmp_N>512) tmp_opti2 = 512
        if(maxval(tmp_int(:))>128) tmp_opti3 = 128

        ! Thread block dimensions
        plan%t_tdma    = dim3(tmp_opti1,1,1)
        plan%t_rdtdma  = dim3(tmp_opti2,1,1)
        plan%t_pack    = dim3(tmp_opti3,1,1)
        plan%t_unpack  = dim3(tmp_opti3,1,1)
        plan%t_intiAb  = dim3(tmp_opti2,1,1)
        
        ! Grid dimensions computed from problem sizes and block sizes
        plan%b_tdma    = dim3(ceiling(dble(Nsys)              /dble(plan%t_tdma  %x)),       1,1)
        plan%b_rdtdma  = dim3(ceiling(dble(tmp_N)             /dble(plan%t_rdtdma%x)),       1,1)
        plan%b_pack    = dim3(ceiling(dble(maxval(tmp_int(:)))/dble(plan%t_pack  %x)),       2,1)
        plan%b_unpack  = dim3(ceiling(dble(maxval(tmp_int(:)))/dble(plan%t_unpack%x)),       2,1)
        plan%b_intiAb  = dim3(ceiling(dble(tmp_N)             /dble(plan%t_intiAb%x)),2*nprocs,1)
       
        deallocate(tmp_int)
    end subroutine pascal_plan_create

    !===================================================================
    ! pascal_plan_clean
    ! -----------------
    ! Release all host/device arrays associated with a given plan.
    ! This should be called once the solver is no longer needed.
    !===================================================================
    subroutine pascal_plan_clean(plan)
        implicit none
        type(ptdma_plan_cuda) :: plan 

        ! Device work arrays
        deallocate(plan%Ard,plan%Brd,plan%Crd,plan%Drd)
        deallocate(plan%Atr,plan%Btr)
        deallocate(plan%Ctr,plan%Dtr)

        ! Host/device gather descriptors
        deallocate(plan%gather_Nrd_local,plan%gather_Ntr_local)
        deallocate(plan%gather_Nrd_local_d,plan%gather_Ntr_local_d)

        deallocate(plan%gather_Nrd_start,plan%gather_Ntr_start)
        deallocate(plan%gather_Nrd_start_d,plan%gather_Ntr_start_d)

        ! Communication buffer descriptors and buffers
        deallocate(plan%bufsubsize_A,plan%bufsubsize_B)
        deallocate(plan%bufstart_A,plan%bufstart_B)
        deallocate(plan%BIGbufsubsize_A,plan%BIGbufsubsize_B)
        deallocate(plan%BIGbufstart_A,plan%BIGbufstart_B)
        deallocate(plan%BIGbuf_A,plan%BIGbuf_B)
    end subroutine pascal_plan_clean

    !===================================================================
    ! pascal_setcudathread
    ! --------------------
    ! (Placeholder)
    ! Intended for externally overriding CUDA block/thread settings
    ! after plan creation. Currently not implemented.
    !===================================================================
    subroutine pascal_setcudathread(plan,int_tdma, int_rdtdma, int_pack, int_unpack, int_intiAb)
        implicit none
        type(ptdma_plan_cuda) :: plan 
        type(dim3) :: int_tdma,int_rdtdma,int_pack,int_unpack,int_intiAb
        
    end subroutine pascal_setcudathread

    !===================================================================
    ! pascal_solver
    ! -------------
    ! Top-level driver for the PaScaL-TDMA algorithm.
    !
    ! For nprocs == 1:
    !   - Directly call tdma_many_cuda on all independent systems.
    !
    ! For nprocs > 1:
    !   1) Apply modified TDMA locally to compress each line
    !   2) Pack reduced coefficients into send buffers (pascalpack)
    !   3) Perform all-to-all communication (pascal_a2av)
    !   4) Unpack to build globally reduced system (pascalunpack)
    !   5) Initialize reduced system RHS (pascalintAb)
    !   6) Solve reduced system using tdma_many_cuda
    !   7) Pack/unpack solutions to distribute interface values
    !   8) Update full solution along each line (pascal_update)
    !
    ! A, B, C, D:
    !   Input/Output device arrays holding tridiagonal coefficients and RHS.
    !===================================================================
    subroutine pascal_solver(plan,A,B,C,D,Nsys,Nrow)
        implicit none
        type(ptdma_plan_cuda) :: plan 
        integer :: Nsys,Nrow
        real*8, device :: A(0:Nsys-1,0:Nrow-1),B(0:Nsys-1,0:Nrow-1),C(0:Nsys-1,0:Nrow-1),D(0:Nsys-1,0:Nrow-1)

        integer :: i
        
        if(plan%nprocs==1) then
            ! Single-process case: standard TDMA on each system
            call tdma_many_cuda<<<plan%b_tdma,plan%t_tdma>>>(A,B,C,D, Nsys, Nrow)
        else
            !----------------------------------------------------------------
            ! Local modified TDMA: compress each line into boundary system
            !----------------------------------------------------------------
            call tdma_modified_cuda<<<plan%b_tdma,plan%t_tdma>>>(A,B,C,D, plan%Ard,plan%Brd,plan%Crd,plan%Drd, Nsys, Nrow)
            
            !----------------------------------------------------------------
            ! PACK: extract reduced coefficients (Ard,Crd,Drd) into BIGbuf_A
            !----------------------------------------------------------------
            do i = 0, plan%nprocs-1          
                call pascalpack<<<plan%b_pack,plan%t_pack>>>(plan%Ard,plan%Nrd_global(0),plan%Nrd_global(1)                  &
                                                            ,plan%gather_Nrd_local_d(0:1,i),plan%gather_Nrd_start_d(0:1,i) &
                                                            ,plan%BIGbuf_A,sum(plan%BIGbufsubsize_A(:))                    &
                                                            ,plan%BIGbufstart_A(i)+0*plan%bufsubsize_A(i)                  )
                call pascalpack<<<plan%b_pack,plan%t_pack>>>(plan%Crd,plan%Nrd_global(0),plan%Nrd_global(1)                  &
                                                            ,plan%gather_Nrd_local_d(0:1,i),plan%gather_Nrd_start_d(0:1,i) &
                                                            ,plan%BIGbuf_A,sum(plan%BIGbufsubsize_A(:))                    &
                                                            ,plan%BIGbufstart_A(i)+1*plan%bufsubsize_A(i)                  )
                call pascalpack<<<plan%b_pack,plan%t_pack>>>(plan%Drd,plan%Nrd_global(0),plan%Nrd_global(1)                  &
                                                            ,plan%gather_Nrd_local_d(0:1,i),plan%gather_Nrd_start_d(0:1,i) &
                                                            ,plan%BIGbuf_A,sum(plan%BIGbufsubsize_A(:))                    &
                                                            ,plan%BIGbufstart_A(i)+2*plan%bufsubsize_A(i)                  )
            end do

            !----------------------------------------------------------------
            ! All-to-all communication of reduced coefficients
            !----------------------------------------------------------------
            call pascal_a2av( plan%BIGbuf_A,sum(plan%BIGbufsubsize_A(:)) &
                            , plan%BIGbufsubsize_A,plan%BIGbufstart_A &
                            , plan%BIGbuf_B,sum(plan%BIGbufsubsize_B(:)) &
                            , plan%BIGbufsubsize_B,plan%BIGbufstart_B &
                            , plan%nprocs, plan%ptdma_world)
    
            !----------------------------------------------------------------
            ! UNPACK: assemble globally reduced system (Atr,Ctr,Dtr)
            !----------------------------------------------------------------
            do i = 0, plan%nprocs-1          
                call pascalunpack<<<plan%b_pack,plan%t_pack>>>(plan%Atr,plan%Ntr_global(0),plan%Ntr_global(1)                  &
                                                            ,plan%gather_Ntr_local_d(0:1,i),plan%gather_Ntr_start_d(0:1,i) &
                                                            ,plan%BIGbuf_B,sum(plan%BIGbufsubsize_B(:))                    &
                                                            ,plan%BIGbufstart_B(i)+0*plan%bufsubsize_B(i)                  )
                call pascalunpack<<<plan%b_pack,plan%t_pack>>>(plan%Ctr,plan%Ntr_global(0),plan%Ntr_global(1)                  &
                                                            ,plan%gather_Ntr_local_d(0:1,i),plan%gather_Ntr_start_d(0:1,i) &
                                                            ,plan%BIGbuf_B,sum(plan%BIGbufsubsize_B(:))                    &
                                                            ,plan%BIGbufstart_B(i)+1*plan%bufsubsize_B(i)                  )
                call pascalunpack<<<plan%b_pack,plan%t_pack>>>(plan%Dtr,plan%Ntr_global(0),plan%Ntr_global(1)                  &
                                                            ,plan%gather_Ntr_local_d(0:1,i),plan%gather_Ntr_start_d(0:1,i) &
                                                            ,plan%BIGbuf_B,sum(plan%BIGbufsubsize_B(:))                    &
                                                            ,plan%BIGbufstart_B(i)+2*plan%bufsubsize_B(i)                  )
            end do

            ! Initialize Btr (reduced RHS) to identity/ones pattern
            call pascalintAb<<<plan%b_intiAb, plan%t_intiAb>>>(plan%Btr, plan%Ntr_global(0), plan%Ntr_global(1))
            
            ! Solve the globally reduced system (size tmp_N x 2*nprocs)
            call tdma_many_cuda<<<plan%b_rdtdma,plan%t_rdtdma>>>(plan%Atr,plan%Btr,plan%Ctr,plan%Dtr &
                                                               ,plan%Ntr_global(0), plan%Ntr_global(1))

            !----------------------------------------------------------------
            ! PACK: gather reduced solutions back into BIGbuf_B
            !----------------------------------------------------------------
            do i = 0, plan%nprocs-1          
                call pascalpack<<<plan%b_pack,plan%t_pack>>>(plan%Dtr,plan%Ntr_global(0),plan%Ntr_global(1)                  &
                                                             ,plan%gather_Ntr_local_d(0:1,i),plan%gather_Ntr_start_d(0:1,i) &
                                                             ,plan%BIGbuf_B,sum(plan%bufsubsize_B(:))                    &
                                                             ,plan%bufstart_B(i)                  )
            end do

            ! All-to-all to distribute interface solutions to all ranks
            call pascal_a2av( plan%BIGbuf_B,sum(plan%bufsubsize_B(:)) &
                            , plan%bufsubsize_B,plan%bufstart_B &
                            , plan%BIGbuf_A,sum(plan%bufsubsize_A(:)) &
                            , plan%bufsubsize_A,plan%bufstart_A &
                            , plan%nprocs, plan%ptdma_world)

            ! UNPACK: scatter reduced solutions into Drd (local reduced RHS)
            do i = 0, plan%nprocs-1          
                call pascalunpack<<<plan%b_pack,plan%t_pack>>>(plan%Drd,plan%Nrd_global(0),plan%Nrd_global(1)                  &
                                                             ,plan%gather_Nrd_local_d(0:1,i),plan%gather_Nrd_start_d(0:1,i) &
                                                             ,plan%BIGbuf_A,sum(plan%BIGbufsubsize_A(:))                    &
                                                             ,plan%bufstart_A(i)                  )
            end do

            ! Final update: reconstruct full solution on each rank
            call pascal_update<<<plan%b_tdma,plan%t_tdma>>>(A,B,C,D, plan%Drd, Nsys, Nrow)
            
        endif
    end subroutine pascal_solver

    !===================================================================
    ! pascal_solver_profiled
    ! ----------------------
    ! Profiled solver entry point for Study drivers. It preserves the same
    ! numerical flow as pascal_solver, but inserts synchronizations around
    ! major CUDA/MPI phases so the measured phase times are meaningful.
    !
    ! This routine is intended for instrumentation and comparison studies.
    ! The regular pascal_solver path is kept free of profiling overhead.
    !===================================================================
    subroutine pascal_solver_profiled(plan,A,B,C,D,Nsys,Nrow,timing)
        implicit none
        type(ptdma_plan_cuda) :: plan
        type(ptdma_timing_cuda) :: timing
        integer :: Nsys,Nrow
        real*8, device :: A(0:Nsys-1,0:Nrow-1),B(0:Nsys-1,0:Nrow-1)
        real*8, device :: C(0:Nsys-1,0:Nrow-1),D(0:Nsys-1,0:Nrow-1)

        integer :: i, ierr
        real(8) :: total_start, phase_start

        call pascal_timing_reset(timing)
        ierr = CudaDeviceSynchronize()
        total_start = MPI_WTIME()

        if(plan%nprocs==1) then
            phase_start = MPI_WTIME()
            call tdma_many_cuda<<<plan%b_tdma,plan%t_tdma>>>(A,B,C,D, Nsys, Nrow)
            ierr = CudaDeviceSynchronize()
            timing%local_compute = MPI_WTIME() - phase_start
            timing%total = MPI_WTIME() - total_start
        else
            phase_start = MPI_WTIME()
            call tdma_modified_cuda<<<plan%b_tdma,plan%t_tdma>>>(A,B,C,D, &
                                                                 plan%Ard,plan%Brd,plan%Crd,plan%Drd, &
                                                                 Nsys, Nrow)
            ierr = CudaDeviceSynchronize()
            timing%local_compute = MPI_WTIME() - phase_start

            phase_start = MPI_WTIME()
            do i = 0, plan%nprocs-1
                call pascalpack<<<plan%b_pack,plan%t_pack>>>(plan%Ard,plan%Nrd_global(0),plan%Nrd_global(1) &
                                                            ,plan%gather_Nrd_local_d(0:1,i)                 &
                                                            ,plan%gather_Nrd_start_d(0:1,i)                 &
                                                            ,plan%BIGbuf_A,sum(plan%BIGbufsubsize_A(:))     &
                                                            ,plan%BIGbufstart_A(i)+0*plan%bufsubsize_A(i)   )
                call pascalpack<<<plan%b_pack,plan%t_pack>>>(plan%Crd,plan%Nrd_global(0),plan%Nrd_global(1) &
                                                            ,plan%gather_Nrd_local_d(0:1,i)                 &
                                                            ,plan%gather_Nrd_start_d(0:1,i)                 &
                                                            ,plan%BIGbuf_A,sum(plan%BIGbufsubsize_A(:))     &
                                                            ,plan%BIGbufstart_A(i)+1*plan%bufsubsize_A(i)   )
                call pascalpack<<<plan%b_pack,plan%t_pack>>>(plan%Drd,plan%Nrd_global(0),plan%Nrd_global(1) &
                                                            ,plan%gather_Nrd_local_d(0:1,i)                 &
                                                            ,plan%gather_Nrd_start_d(0:1,i)                 &
                                                            ,plan%BIGbuf_A,sum(plan%BIGbufsubsize_A(:))     &
                                                            ,plan%BIGbufstart_A(i)+2*plan%bufsubsize_A(i)   )
            end do
            ierr = CudaDeviceSynchronize()
            timing%pack_forward = MPI_WTIME() - phase_start

            phase_start = MPI_WTIME()
            call pascal_a2av( plan%BIGbuf_A,sum(plan%BIGbufsubsize_A(:)),plan%BIGbufsubsize_A,plan%BIGbufstart_A &
                            , plan%BIGbuf_B,sum(plan%BIGbufsubsize_B(:)),plan%BIGbufsubsize_B,plan%BIGbufstart_B &
                            , plan%nprocs, plan%ptdma_world)
            timing%mpi_forward = MPI_WTIME() - phase_start

            phase_start = MPI_WTIME()
            do i = 0, plan%nprocs-1
                call pascalunpack<<<plan%b_pack,plan%t_pack>>>(plan%Atr,plan%Ntr_global(0),plan%Ntr_global(1) &
                                                              ,plan%gather_Ntr_local_d(0:1,i)                 &
                                                              ,plan%gather_Ntr_start_d(0:1,i)                 &
                                                              ,plan%BIGbuf_B,sum(plan%BIGbufsubsize_B(:))     &
                                                              ,plan%BIGbufstart_B(i)+0*plan%bufsubsize_B(i)   )
                call pascalunpack<<<plan%b_pack,plan%t_pack>>>(plan%Ctr,plan%Ntr_global(0),plan%Ntr_global(1) &
                                                              ,plan%gather_Ntr_local_d(0:1,i)                 &
                                                              ,plan%gather_Ntr_start_d(0:1,i)                 &
                                                              ,plan%BIGbuf_B,sum(plan%BIGbufsubsize_B(:))     &
                                                              ,plan%BIGbufstart_B(i)+1*plan%bufsubsize_B(i)   )
                call pascalunpack<<<plan%b_pack,plan%t_pack>>>(plan%Dtr,plan%Ntr_global(0),plan%Ntr_global(1) &
                                                              ,plan%gather_Ntr_local_d(0:1,i)                 &
                                                              ,plan%gather_Ntr_start_d(0:1,i)                 &
                                                              ,plan%BIGbuf_B,sum(plan%BIGbufsubsize_B(:))     &
                                                              ,plan%BIGbufstart_B(i)+2*plan%bufsubsize_B(i)   )
            end do
            ierr = CudaDeviceSynchronize()
            timing%unpack_forward = MPI_WTIME() - phase_start

            phase_start = MPI_WTIME()
            call pascalintAb<<<plan%b_intiAb, plan%t_intiAb>>>(plan%Btr, plan%Ntr_global(0), plan%Ntr_global(1))
            call tdma_many_cuda<<<plan%b_rdtdma,plan%t_rdtdma>>>(plan%Atr,plan%Btr,plan%Ctr,plan%Dtr, &
                                                                 plan%Ntr_global(0), plan%Ntr_global(1))
            ierr = CudaDeviceSynchronize()
            timing%reduced_compute = MPI_WTIME() - phase_start

            phase_start = MPI_WTIME()
            do i = 0, plan%nprocs-1
                call pascalpack<<<plan%b_pack,plan%t_pack>>>(plan%Dtr,plan%Ntr_global(0),plan%Ntr_global(1) &
                                                            ,plan%gather_Ntr_local_d(0:1,i)                 &
                                                            ,plan%gather_Ntr_start_d(0:1,i)                 &
                                                            ,plan%BIGbuf_B,sum(plan%bufsubsize_B(:))        &
                                                            ,plan%bufstart_B(i)                             )
            end do
            ierr = CudaDeviceSynchronize()
            timing%pack_backward = MPI_WTIME() - phase_start

            phase_start = MPI_WTIME()
            call pascal_a2av( plan%BIGbuf_B,sum(plan%bufsubsize_B(:)),plan%bufsubsize_B,plan%bufstart_B &
                            , plan%BIGbuf_A,sum(plan%bufsubsize_A(:)),plan%bufsubsize_A,plan%bufstart_A &
                            , plan%nprocs, plan%ptdma_world)
            timing%mpi_backward = MPI_WTIME() - phase_start

            phase_start = MPI_WTIME()
            do i = 0, plan%nprocs-1
                call pascalunpack<<<plan%b_pack,plan%t_pack>>>(plan%Drd,plan%Nrd_global(0),plan%Nrd_global(1) &
                                                              ,plan%gather_Nrd_local_d(0:1,i)                 &
                                                              ,plan%gather_Nrd_start_d(0:1,i)                 &
                                                              ,plan%BIGbuf_A,sum(plan%BIGbufsubsize_A(:))     &
                                                              ,plan%bufstart_A(i)                             )
            end do
            ierr = CudaDeviceSynchronize()
            timing%unpack_backward = MPI_WTIME() - phase_start

            phase_start = MPI_WTIME()
            call pascal_update<<<plan%b_tdma,plan%t_tdma>>>(A,B,C,D, plan%Drd, Nsys, Nrow)
            ierr = CudaDeviceSynchronize()
            timing%update_compute = MPI_WTIME() - phase_start

            timing%total = MPI_WTIME() - total_start
        endif
    end subroutine pascal_solver_profiled

    !===================================================================
    ! para
    ! ----
    ! Simple 1D block partition of a global index range [nsta, nend]
    ! into "nprocs" contiguous subranges. This rank (myrank) receives
    ! the subrange [indx_a, indx_b].
    !===================================================================
    subroutine para(nsta,nend,nprocs,myrank,indx_a,indx_b)
        implicit none
        integer :: nsta,nend,nprocs,myrank,indx_a,indx_b
        integer :: iwork1, iwork2
        
        iwork1 = int((nend-nsta+1)/nprocs)
        iwork2 = mod((nend-nsta+1),nprocs)
        indx_a = myrank*iwork1 + nsta +min(myrank,iwork2)
        indx_b = indx_a + iwork1 -1
        if(iwork2 > myrank) indx_b = indx_b +1
        
    end subroutine para

    !===================================================================
    ! pascalpack
    ! ----------
    ! CUDA kernel that packs a 2D sub-block of A into a linear buffer.
    !   - A        : input 2D array (device)
    !   - pack_*   : size and starting indices of the sub-block
    !   - buf_A    : output 1D buffer on device
    !   - bufpoint : offset inside buf_A where this sub-block is stored
    ! Used to assemble send buffers for MPI all-to-all communication.
    !===================================================================
    attributes(global) subroutine pascalpack(A,n1,n2,pack_subsize,pack_start,buf_A,bufsize,bufpoint)
        use cudafor
        implicit none
        integer,value :: n1,n2,bufsize,bufpoint
        integer, device :: pack_subsize(0:1),pack_start(0:1)
        real*8, device :: A(0:n1-1,0:n2-1)
        real*8, device :: buf_A(0:bufsize-1)
        
        integer :: i,j,k,ierr
        integer :: indexi,indexj,indexk,indexbf
        
        ! Compute (i,j) index from thread / block indices
        i = (blockidx%x - 1)*blockdim%x + (threadidx%x-1)
        j = (blockidx%y - 1)*blockdim%y + (threadidx%y-1)

        ! Global (row, column) indices in A for this thread
        indexi = i + pack_start(0)
        indexj = j + pack_start(1)
        
        ! If within bounds, write into contiguous buffer
        if ( (indexi<n1)   &
        .and.(indexj<n2)  ) then
            indexbf =i + j*pack_subsize(0) + bufpoint
            buf_A(indexbf) = A(indexi,indexj)
        end if        
        
    end subroutine pascalpack
    
    !===================================================================
    ! pascalunpack
    ! ------------
    ! CUDA kernel that unpacks data from a 1D buffer back into a 2D
    ! sub-block of array A. This is the inverse of pascalpack.
    !===================================================================
    attributes(global) subroutine pascalunpack(A,n1,n2,pack_subsize,pack_start,buf_A,bufsize,bufpoint)
        use cudafor
        implicit none
        integer,value :: n1,n2,bufsize,bufpoint
        integer, device :: pack_subsize(0:1),pack_start(0:1)
        real*8, device :: A(0:n1-1,0:n2-1)
        real*8, device :: buf_A(0:bufsize-1)
        
        integer :: i,j,k,ierr
        integer :: indexi,indexj,indexk,indexbf
        
        ! Compute (i,j) index from thread / block indices
        i = (blockidx%x - 1)*blockdim%x + (threadidx%x-1)
        j = (blockidx%y - 1)*blockdim%y + (threadidx%y-1)

        ! Global (row, column) indices in A for this thread
        indexi = i + pack_start(0)
        indexj = j + pack_start(1)
        
        ! If within bounds, read from contiguous buffer into A
        if ( (indexi<n1)   &
        .and.(indexj<n2)  ) then
            indexbf =i + j*pack_subsize(0) + bufpoint
            A(indexi,indexj) = buf_A(indexbf) 
        end if        
        
    end subroutine pascalunpack

    
    !===================================================================
    ! tdma_many_cuda
    ! --------------
    ! CUDA kernel implementing the standard Thomas algorithm for many
    ! independent tridiagonal systems in parallel.
    !
    !   - Each thread i solves one system (over j = 0..nrow-1)
    !   - a,b,c : tridiagonal coefficients
    !   - d     : RHS; overwritten in-place by the solution
    !===================================================================
    attributes(global) subroutine tdma_many_cuda(a,b,c,d, nsys, nrow)
        integer, value :: nsys, nrow
        real*8, device :: a(0:nsys-1,0:nrow-1),b(0:nsys-1,0:nrow-1),c(0:nsys-1,0:nrow-1),d(0:nsys-1,0:nrow-1)

        real*8:: a1_sh
        real*8:: b1_sh
        real*8:: c1_sh,c0_sh
        real*8:: d1_sh,d0_sh

        integer :: i,j
        integer :: ti
        real*8  :: r

        ti = (threadidx%x-1)
        i  = (threadidx%x-1) + (blockidx%x-1)*blockdim%x
    
        if(i<nsys) then
            ! Forward sweep: eliminate lower diagonal
            b1_sh = b(i,0)
            c1_sh = c(i,0)
            d1_sh = d(i,0)

            d1_sh = d1_sh/b1_sh
            c1_sh = c1_sh/b1_sh

            d(i,0)=d1_sh
            c(i,0)=c1_sh

            do j=1, nrow-1
                c0_sh = c1_sh
                d0_sh = d1_sh

                a1_sh = a(i,j) 
                b1_sh = b(i,j) 
                c1_sh = c(i,j) 
                d1_sh = d(i,j)
                
                r = 1.d0/(b1_sh-a1_sh*c0_sh)
                d1_sh = r*(d1_sh-a1_sh*d0_sh)
                c1_sh = r*c1_sh

                d(i,j) = d1_sh
                c(i,j) = c1_sh
            enddo

            ! Back substitution
            do j=nrow-2,0,-1
                c0_sh = c(i,j)
                d0_sh = d(i,j)
                d0_sh = d0_sh-c0_sh*d1_sh
                d1_sh = d0_sh
                d(i,j) = d0_sh
            enddo
        endif
    end subroutine tdma_many_cuda

    !===================================================================
    ! tdma_modified_cuda
    ! ------------------
    ! Variant of TDMA that:
    !   - Applies a local transformation to compress each system into
    !     a 2×2 reduced representation (stored in a_rd,b_rd,c_rd,d_rd)
    !   - These reduced systems are later coupled across ranks and
    !     solved as a global reduced problem.
    !===================================================================
    attributes(global) subroutine tdma_modified_cuda(a,b,c,d, a_rd,b_rd,c_rd,d_rd, nsys, nrow)
        integer, value :: nsys, nrow
        real*8, device :: a(0:nsys-1,0:nrow-1),b(0:nsys-1,0:nrow-1),c(0:nsys-1,0:nrow-1),d(0:nsys-1,0:nrow-1)
        real*8, device :: a_rd(0:nsys-1,0:1),b_rd(0:nsys-1,0:1),c_rd(0:nsys-1,0:1),d_rd(0:nsys-1,0:1)

        real*8:: a1_sh,a0_sh
        real*8:: b1_sh,b0_sh
        real*8:: c1_sh,c0_sh
        real*8:: d1_sh,d0_sh

        integer :: i,j
        integer :: ti
        real*8  :: r,r0_sh

        ti = (threadidx%x-1)
        i  = (threadidx%x-1) + (blockidx%x-1)*blockdim%x
    
        if(i<nsys) then
            ! First two rows: initial normalization
            a0_sh = a(i,0)
            b0_sh = b(i,0)
            c0_sh = c(i,0)
            d0_sh = d(i,0)

            a0_sh = a0_sh/b0_sh
            c0_sh = c0_sh/b0_sh
            d0_sh = d0_sh/b0_sh

            a(i,0)=a0_sh
            c(i,0)=c0_sh
            d(i,0)=d0_sh
            
            a1_sh = a(i,1)
            b1_sh = b(i,1)
            c1_sh = c(i,1)
            d1_sh = d(i,1)

            a1_sh = a1_sh/b1_sh
            c1_sh = c1_sh/b1_sh
            d1_sh = d1_sh/b1_sh

            a(i,1)=a1_sh
            c(i,1)=c1_sh
            d(i,1)=d1_sh

            ! Forward elimination for rows 2..nrow-1 with modified scheme
            do j=2, nrow-1
                a0_sh = a1_sh
                c0_sh = c1_sh
                d0_sh = d1_sh

                a1_sh = a(i,j) 
                b1_sh = b(i,j) 
                c1_sh = c(i,j) 
                d1_sh = d(i,j)
                
                r0_sh = 1.d0/(b1_sh-a1_sh*c0_sh)
                d1_sh = r0_sh*(d1_sh-a1_sh*d0_sh)
                c1_sh = r0_sh*c1_sh
                a1_sh =-r0_sh*a1_sh*a0_sh

                a(i,j) = a1_sh
                c(i,j) = c1_sh
                d(i,j) = d1_sh
            enddo

            ! Store reduced representation (end row)
            a_rd(i,1) = a1_sh
            b_rd(i,1) = 1.d0
            c_rd(i,1) = c1_sh
            d_rd(i,1) = d1_sh

            ! Backward sweep to compress toward the first row
            a1_sh = a0_sh
            c1_sh = c0_sh
            d1_sh = d0_sh

            do j=nrow-3,1,-1
                a0_sh = a(i,j)
                d0_sh = d(i,j)
                c0_sh = c(i,j)

                a0_sh = a0_sh-c0_sh*a1_sh
                d0_sh = d0_sh-c0_sh*d1_sh
                c0_sh =-c0_sh*c1_sh
                
                a1_sh = a0_sh
                d1_sh = d0_sh
                c1_sh = c0_sh

                a(i,j) = a0_sh
                d(i,j) = d0_sh
                c(i,j) = c0_sh
            enddo

            ! Final reduction to a 2×2 system stored in a_rd, b_rd, c_rd, d_rd
            a0_sh = a(i,0)
            d0_sh = d(i,0)
            c0_sh = c(i,0)

            r0_sh = 1.d0/(1.d0-a1_sh*c0_sh)
            a0_sh = r0_sh*a0_sh
            d0_sh = r0_sh*(d0_sh-c0_sh*d1_sh)
            c0_sh =-r0_sh*c0_sh*c1_sh

            a(i,0) = a0_sh
            d(i,0) = d0_sh
            c(i,0) = c0_sh

            a_rd(i,0) = a0_sh
            b_rd(i,0) = 1.d0
            c_rd(i,0) = c0_sh
            d_rd(i,0) = d0_sh
        endif
    end subroutine tdma_modified_cuda

    !===================================================================
    ! pascal_update
    ! -------------
    ! CUDA kernel that reconstructs the full solution along each line
    ! using the reduced solution at both ends (d_rd):
    !   - ds_sh = solution at start
    !   - de_sh = solution at end
    !   - Interior points are adjusted by subtracting contributions
    !     from a(i,j)*ds_sh and c(i,j)*de_sh.
    !===================================================================
    attributes(global) subroutine pascal_update(a,b,c,d, d_rd, nsys, nrow)
        integer, value :: nsys, nrow
        real*8, device :: a(0:nsys-1,0:nrow-1),b(0:nsys-1,0:nrow-1),c(0:nsys-1,0:nrow-1),d(0:nsys-1,0:nrow-1)
        real*8, device :: d_rd(0:nsys-1,0:1)

        real*8:: ds_sh,de_sh

        integer :: i,j
        integer :: ti
        real*8  :: r,r0_sh

        ti = (threadidx%x-1)
        i  = (threadidx%x-1) + (blockidx%x-1)*blockdim%x
        if(i<nsys) then
            ds_sh =d_rd(i,0)
            de_sh =d_rd(i,1)
            !call synthreads()

            ! Apply boundary solutions at both ends
            d(i,0)      = ds_sh
            d(i,nrow-1) = de_sh

            ! Correct interior entries using boundary contributions
            do j=1,nrow-2
                d(i,j) = d(i,j) - a(i,j)*ds_sh - c(i,j)*de_sh
            enddo
        endif
    end subroutine pascal_update
    
    !===================================================================
    ! pascal_a2av
    ! -----------
    ! Wrapper around MPI_ALLTOALLV for double precision buffers.
    ! A and B are device arrays; MPI is expected to be CUDA-aware
    ! so that device pointers can be passed directly.
    !
    !   A        : send buffer (device)
    !   B        : receive buffer (device)
    !   send*    : counts/displacements for A
    !   recv*    : counts/displacements for B
    !===================================================================
    subroutine pascal_a2av(A,Asize,sendcount,senddisp, B,Bsize,recvcount,recvdisp, nprocs,communicator)
        use mpi
        use cudafor
        implicit none
        integer:: Asize,Bsize
        integer :: nprocs
        real*8, dimension(:), device :: A(0:Asize-1),B(0:Bsize-1) 
        integer, dimension(:) :: sendcount(0:nprocs-1),senddisp(0:nprocs-1)
        integer, dimension(:) :: recvcount(0:nprocs-1),recvdisp(0:nprocs-1)
        integer :: communicator

        integer :: i,myrank,ierr   
        integer :: requestA(0:1024-1),requestB(0:1024-1)     

        ierr = CudaDeviceSynchronize()
        call MPI_ALLTOALLV(A,sendcount,senddisp,MPI_DOUBLE, B,recvcount,recvdisp,MPI_DOUBLE, communicator, ierr)
        ! call MPI_COMM_RANK(communicator, myrank, ierr)
        ! do i = 0, nprocs-1
        !     call MPI_ISEND(A(senddisp(i)),sendcount(i),MPI_DOUBLE,i,111,communicator,requestA(i),ierr)
        ! end do
        
        ! do i = 0, nprocs-1
        !     call MPI_IRECV(B(recvdisp(i)),recvcount(i),MPI_DOUBLE,i,111,communicator,requestB(i),ierr)
        ! end do
        
        ! call MPI_WAITALL(nprocs, requestA, MPI_STATUSES_IGNORE, ierr)
        ! call MPI_WAITALL(nprocs, requestB, MPI_STATUSES_IGNORE, ierr)
        
    end subroutine pascal_a2av

    !===================================================================
    ! pascalintAb
    ! -----------
    ! Simple CUDA kernel that initializes a 2D array A to 1.0 inside
    ! the given bounds. Used to set up a reduced system RHS or
    ! identity-like structure in the transformed space.
    !===================================================================
    attributes(global) subroutine pascalintAb(A,in_n,in_m)
        use cudafor
        implicit none
        integer, value :: in_n,in_m
        real*8, device :: A(0:in_n-1,0:in_m-1)

        integer :: i,j,k,ierr
        
        ! Compute (i,j) from block/thread indices
        i = (blockidx%x - 1)*blockdim%x + (threadidx%x-1)
        j = (blockidx%y - 1)*blockdim%y + (threadidx%y-1)
        
        ! Initialize within bounds
        if ( (i<in_n) .and.(j<in_m)   ) then
            A(i,j) = 1.d0
        end if        
        
    end subroutine pascalintAb

end module PaScaL_TDMA_cuda
