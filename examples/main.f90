!!======================================================================================================================
!> @file        main.f90
!> @brief       This file contains the main subroutines for an example problem of PaScaL_TDMA.
!> @details     The target example problem is a three-dimensional(3D) time-dependent heat conduction problem 
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
!> @brief       Main execution program for the example problem of PaScaL_TDMA.
!>
program main

    use mpi
    use global
    use mpi_subdomain
    use mpi_topology
    use solve_theta

#ifdef _CUDA
    use cudafor
    use solve_theta, only : solve_theta_plan_many_cuda
    use nvtx
#else
    use PaScaL_TDMA
#endif
    
    implicit none
 
    integer :: nprocs, myrank   ! Number of MPI processes and rank ID in MPI_COMM_WORLD
    integer :: ierr
    double precision, allocatable, dimension(:, :, :) :: theta_sub  ! Main 3-D variable to be solved

#ifdef _CUDA
    integer :: local_rank, pvd, istat, nDevices
#endif

    call MPI_Init(ierr)
    call MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr)
    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr)
    
#ifdef _CUDA
    istat = cudaGetDeviceCount(nDevices) ! Count only GPUs in single node
    local_rank = mod(myrank,nDevices) ! and devide it for GPU affinity.
    istat = cudaSetDevice(local_rank)
#endif

    if(myrank==0) write(*,*) '[Main] The main simulation starts! '
    ! Periodicity in the simulation domain
    period(0)=.true.; period(1)=.false.; period(2)=.true.
    
    ! Read PARA_INPUT.inp file and setup the parameters for the simulation domain and conditions.
    call global_inputpara(np_dim)

    ! Create cartesian MPI topology and build subcommunicator in the cartesian topology.
    call mpi_topology_make()

    ! Divide the simulation domain into sub-domains.
    call mpi_subdomain_make(comm_1d_x%nprocs, comm_1d_x%myrank, &
                                   comm_1d_y%nprocs, comm_1d_y%myrank, &
                                   comm_1d_z%nprocs, comm_1d_z%myrank )

    ! Make ghostcells as communication buffer using derived datatypes.
    call mpi_subdomain_make_ghostcell_ddtype()

    if(myrank==0) write(*,*) '[Main] Subdomains and topology created! '

    ! Allocate the main variable in the sub-domain.
    allocate(theta_sub(0:nx_sub, 0:ny_sub, 0:nz_sub))

    ! Setup the parameters in the sub-domain.
    call mpi_subdomain_indices(comm_1d_y%myrank, comm_1d_y%nprocs)            
    call mpi_subdomain_mesh(comm_1d_x%myrank,comm_1d_y%myrank,comm_1d_z%myrank, &
                                   comm_1d_x%nprocs,comm_1d_y%nprocs,comm_1d_z%nprocs)

    ! Initialize values in the simulation domain and boundary conditions.
    call mpi_subdomain_initialization(theta_sub,comm_1d_y%myrank, comm_1d_y%nprocs)
    call mpi_subdomain_ghostcell_update(theta_sub, comm_1d_x, comm_1d_y, comm_1d_z)
    call mpi_subdomain_boundary(theta_sub, comm_1d_y%myrank, comm_1d_y%nprocs)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    if(myrank==0) write(*,*) '[Main] The preparation for simulation finished! '

    if(myrank==0) write(*,*) '[Main] Solving the 3D heat equation! '

#ifdef _CUDA
    call solve_theta_plan_many_cuda(theta_sub)
#else
    call solve_theta_plan_many(theta_sub)
#endif
 
    if(myrank==0) write(*,*) '[Main] Solving the 3D heat equation complete! '

    ! Write the values of the field variable in an output file if required.
    call field_file_write(myrank, nprocs, theta_sub)
    if(myrank==0) write(*,*) '[Main] Solutions have been written to files! '
   
    ! Deallocate the variables and clean the module variables.
    deallocate(theta_sub)
    call mpi_subdomain_clean()
    call mpi_topology_clean()
    if(myrank==0) write(*,*) '[Main] The main simulation complete! '

    call MPI_Finalize(ierr)
end program main

!>
!> @brief       Write the values of the field variable in an output file.
!> @details     It writes the values of the field variables using two methods; single task IO and post assembly IO.
!>              In the single task IO, each process writes the values in its own file. In the post assembly IO,  
!>              each process send data to a master process, which writes the obtained data in a single file.
!> @param       myrank      Rank ID in MPI_COMM_WORLD
!> @param       nprocs      Number of MPI process in MPI_COMM_WORLD
!> @param       theta_sub   Field variable to write
!>
subroutine field_file_write(myrank, nprocs, theta_sub)

    use mpi
    use global
    use mpi_subdomain
    use mpi_topology

    implicit none

    integer, intent(in) :: myrank, nprocs
    double precision, dimension(0:nx_sub, 0:ny_sub, 0:nz_sub), intent(in) :: theta_sub

    !---- MPI write ----
    character(len=256) :: filename                                  ! Output file name for each process
    double precision, allocatable, dimension(:, :, :) :: theta_all  ! Global field variable to write in a single file 
    integer, allocatable, dimension(:,:) :: cart_coord              ! MPI process coordinates in cartesin topology
    integer :: sizes(0:2), subsizes(0:2), starts(0:2)                           ! For MPI_Type_create_subarray
    integer, allocatable, dimension(:) :: nxm_sub_cnt,  nym_sub_cnt,  nzm_sub_cnt  ! For MPI_Type_create_subarray
    integer, allocatable, dimension(:) :: nxm_sub_disp, nym_sub_disp, nzm_sub_disp ! For MPI_Type_create_subarray
    integer :: ddtype_data_write_send, request_send                             ! DDT from MPI_Type_create_subarray
    integer, allocatable, dimension(:) :: ddtype_data_write_recv, request_recv  ! DDT from MPI_Type_create_subarray
    integer :: i, j, k, ierr
    integer :: status(MPI_STATUS_SIZE)

    ! Singe task IO: each process writes to its own output file.
    write (filename, '("mpi_Tfield_sub_af",5I1,".PLT" )' ) int(myrank/10000)-int(myrank/100000)*10 &
                                                                &,int(myrank/1000 )-int(myrank/10000 )*10 &
                                                                &,int(myrank/100  )-int(myrank/1000  )*10 &
                                                                &,int(myrank/10   )-int(myrank/100   )*10 &
                                                                &,int(myrank/1    )-int(myrank/10    )*10 
    open(unit=myrank,file=filename)
    
    write(myrank,*) 'VARIABLES="X","Y","Z","THETA"'
    write(myrank,*) 'zone t="',1,'"','i=',nx_sub+1,'j=',ny_sub+1,'k=',nz_sub+1

    do k=1,nz_sub-1
        do j=1,ny_sub-1
            do i=1,nx_sub-1
                write(myrank,'(3D14.6,1E22.14)') X_sub(i), Y_sub(j), Z_sub(k), theta_sub(i,j,k)
            enddo
        enddo
    enddo

    close(myrank)

    ! Single file IO - process 0 gathers all results from other variables by using derived datatypes.
    allocate(nxm_sub_cnt(0:comm_1d_x%nprocs-1), nxm_sub_disp(0:comm_1d_x%nprocs-1))
    allocate(nym_sub_cnt(0:comm_1d_y%nprocs-1), nym_sub_disp(0:comm_1d_y%nprocs-1))
    allocate(nzm_sub_cnt(0:comm_1d_z%nprocs-1), nzm_sub_disp(0:comm_1d_z%nprocs-1))
    allocate(ddtype_data_write_recv(0:nprocs-1))
    allocate(request_recv(0:nprocs-1))
    allocate(theta_all(1:nxm,1:nym,1:nzm))
    theta_all(:,:,:) = 0.0d0

    ! Building derived datatype for single file IO using the post assembly IO.
    ! Derived datatype for sending data, excluding ghostcells.
    sizes    = (/nx_sub+1,ny_sub+1,nz_sub+1/)
    subsizes = (/nx_sub-1,ny_sub-1,nz_sub-1/)
    starts   = (/1,      1,      1/)

    call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN,  &
                                 MPI_DOUBLE_PRECISION, ddtype_data_write_send, ierr)
    call MPI_Type_commit( ddtype_data_write_send, ierr)

    ! Derived datatype for receiving data and assigning the local data array in the global array for all MPI processes.
    allocate( cart_coord(0:2,0:nprocs-1) )
    do i = 0, nprocs-1
        call MPI_Cart_coords(mpi_world_cart, i, 3, cart_coord(:,i), ierr )
    enddo

    call MPI_Allgather(nx_sub, 1, MPI_INTEGER, nxm_sub_cnt,1, MPI_INTEGER, comm_1d_x%mpi_comm, ierr)
    nxm_sub_cnt(:) = nxm_sub_cnt(:) - 1
    call MPI_Allgather(ny_sub, 1, MPI_INTEGER, nym_sub_cnt,1, MPI_INTEGER, comm_1d_y%mpi_comm, ierr)
    nym_sub_cnt(:) = nym_sub_cnt(:) - 1
    call MPI_Allgather(nz_sub, 1, MPI_INTEGER, nzm_sub_cnt,1, MPI_INTEGER, comm_1d_z%mpi_comm, ierr)
    nzm_sub_cnt(:) = nzm_sub_cnt(:) - 1

    nxm_sub_disp(0) = 0
    do i = 1, comm_1d_x%nprocs-1
        nxm_sub_disp(i)=sum(nxm_sub_cnt(0:i-1))
    enddo

    nym_sub_disp(0) = 0
    do i = 1, comm_1d_y%nprocs-1
        nym_sub_disp(i)=sum(nym_sub_cnt(0:i-1))
    enddo

    nzm_sub_disp(0) = 0
    do i = 1, comm_1d_z%nprocs-1
        nzm_sub_disp(i)=sum(nzm_sub_cnt(0:i-1))
    enddo

    ! Derived datatypes for every process are required for receiving data.
    do i = 0, nprocs-1
        sizes    = (/nxm,    nym,    nzm/)
        subsizes = (/nxm_sub_cnt(cart_coord(0,i)), nym_sub_cnt(cart_coord(1,i)), nzm_sub_cnt(cart_coord(2,i))/)
        starts   = (/nxm_sub_disp(cart_coord(0,i)), nym_sub_disp(cart_coord(1,i)), nzm_sub_disp(cart_coord(2,i))/)

        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                      MPI_DOUBLE_PRECISION, ddtype_data_write_recv(i), ierr)
        call MPI_Type_commit( ddtype_data_write_recv(i), ierr)
    enddo

    ! Each rank sends data to process 0 using DDT.
    call MPI_Isend(theta_sub(0,0,0), 1, ddtype_data_write_send, 0, 101, MPI_COMM_WORLD, request_send, ierr)

    ! Process 0 receives data from all ranks.
    if(myrank == 0 ) then
        do i = 0, nprocs-1
            call MPI_Irecv(theta_all(1,1,1), 1, ddtype_data_write_recv(i), i, 101, MPI_COMM_WORLD, request_recv(i), ierr)
        enddo
    endif

    call MPI_Wait(request_send, status, ierr)
    if(myrank == 0 ) then
        call MPI_Waitall(nprocs, request_recv, MPI_STATUSES_IGNORE, ierr)
    endif

    ! Write gathered data to a single file.
    if(myrank == 0) then
        open(unit=myrank,file="T_field_all.dat")
        do k = 1, nzm
            do j = 1, nym
                do i = 1, nxm
                    write(myrank,'(3D14.6,1E22.14)') dx*dble(i-1),dy*dble(j-1),dz*dble(k-1),theta_all(i,j,k)
                enddo
            enddo
        enddo
        close(myrank)
    endif

    ! Process complete: free the datatypes and deallocate variables.
    call MPI_Type_free(ddtype_data_write_send, ierr)
    do i = 0, nprocs-1
        call MPI_Type_free(ddtype_data_write_recv(i), ierr)
    enddo
    deallocate(theta_all)
    deallocate(nxm_sub_cnt, nxm_sub_disp)
    deallocate(nym_sub_cnt, nym_sub_disp)
    deallocate(nzm_sub_cnt, nzm_sub_disp)
    deallocate(ddtype_data_write_recv)
    deallocate( cart_coord )

end subroutine field_file_write
