program example_fortran_profile
    use mpi
    use cudafor
    use PaScaL_TDMA_cuda
    implicit none

    integer :: ierr, myrank, nprocs
    integer :: ngpu, gpurank
    integer :: n1, n2, n3, iterations
    integer :: nthread_modithomas, nthread_reduced
    integer :: n1sub, n2sub, n3sub, nsys
    integer :: ia, ib, iter
    integer :: nrow_min, nrow_max
    integer, parameter :: n_timing_fields = 13
    real(8) :: phase_local(n_timing_fields)
    real(8) :: phase_max(n_timing_fields), phase_sum(n_timing_fields)

    real(8), allocatable, dimension(:,:,:) :: Aa, Ab, Ac, B
    real(8), allocatable, dimension(:,:,:), device :: Aa_d, Ab_d, Ac_d, B_d
    type(ptdma_plan_cuda) :: plan
    type(ptdma_timing_cuda) :: timing

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

    call parse_positive_arg(1, 64, n1, "n1", myrank)
    call parse_positive_arg(2, 64, n2, "n2", myrank)
    call parse_positive_arg(3, 2048, n3, "n3", myrank)
    call parse_positive_arg(4, 10, iterations, "iterations", myrank)
    call parse_positive_arg(5, 128, nthread_modithomas, "tdma_threads", myrank)
    call parse_positive_arg(6, 128, nthread_reduced, "reduced_threads", myrank)

    if (command_argument_count() > 6) then
        if (myrank == 0) then
            write(*,*) "usage: example_fortran_profile [n1] [n2] [n3]"
            write(*,*) "       [iterations] [tdma_threads] [reduced_threads]"
        endif
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    endif

    ierr = cudaGetDeviceCount(ngpu)
    if (ngpu <= 0) then
        if (myrank == 0) write(*,*) "No CUDA device is visible."
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    endif
    gpurank = mod(myrank, ngpu)
    ierr = cudaSetDevice(gpurank)
    ierr = cudaDeviceSynchronize()

    call para(0, n3 - 1, nprocs, myrank, ia, ib)
    n1sub = n1
    n2sub = n2
    n3sub = ib - ia + 1
    nsys = n1sub * n2sub

    call MPI_REDUCE(n3sub, nrow_min, 1, MPI_INTEGER, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(n3sub, nrow_max, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

    allocate(Aa(0:n1sub-1,0:n2sub-1,0:n3sub-1), Aa_d(0:n1sub-1,0:n2sub-1,0:n3sub-1))
    allocate(Ab(0:n1sub-1,0:n2sub-1,0:n3sub-1), Ab_d(0:n1sub-1,0:n2sub-1,0:n3sub-1))
    allocate(Ac(0:n1sub-1,0:n2sub-1,0:n3sub-1), Ac_d(0:n1sub-1,0:n2sub-1,0:n3sub-1))
    allocate(B (0:n1sub-1,0:n2sub-1,0:n3sub-1), B_d (0:n1sub-1,0:n2sub-1,0:n3sub-1))

    Aa(:,:,:) = 1.0d0
    Ab(:,:,:) = -2.0d0
    Ac(:,:,:) = 1.0d0
    B (:,:,:) = 0.0d0
    if (myrank == 0) B(:,:,0) = -1.0d0
    if (myrank == nprocs - 1) B(:,:,n3sub-1) = -1.0d0

    Aa_d = Aa
    Ab_d = Ab
    Ac_d = Ac
    B_d  = B

    call pascal_plan_create(plan, nsys, MPI_COMM_WORLD, myrank, nprocs, &
                            nthread_modithomas, nthread_reduced)

    if (myrank == 0) then
        write(*,'(A)',advance='no') "solver,implementation,nranks,n1,n2,n3,nsys,"
        write(*,'(A)',advance='no') "nrow_min,nrow_max,iter,iterations,mpi_mode,total_s_max,total_s_avg,"
        write(*,'(A)',advance='no') "local_compute_s_max,pack_forward_s_max,mpi_forward_s_max,"
        write(*,'(A)',advance='no') "unpack_forward_s_max,reduced_compute_s_max,pack_backward_s_max,"
        write(*,'(A)',advance='no') "mpi_backward_s_max,unpack_backward_s_max,update_compute_s_max,"
        write(*,'(A)') "compute_s_max,communication_s_max,packing_s_max"
    endif

    do iter = 0, iterations - 1
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        ! The profiling driver intentionally repeats only the solve call.
        ! Coefficients are not reinitialized between iterations.
        call pascal_solver_profiled(plan, Aa_d, Ab_d, Ac_d, B_d, nsys, n3sub, timing)

        if (iter == 0) then
            B = B_d
            call write_correctness_csv_if_requested(B, n1, n2, n3, nsys, n3sub, &
                                                    nrow_min, nrow_max, ia, ib, &
                                                    nprocs, myrank)
        endif

        call fill_timing_fields(timing, phase_local)
        call MPI_REDUCE(phase_local, phase_max, n_timing_fields, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(phase_local, phase_sum, n_timing_fields, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

        if (myrank == 0) then
            call write_csv_row(nprocs, n1, n2, n3, nsys, nrow_min, nrow_max, &
                               iter, iterations, phase_max, phase_sum(1) / dble(nprocs))
        endif
    enddo

    call pascal_plan_clean(plan)
    deallocate(Aa, Aa_d)
    deallocate(Ab, Ab_d)
    deallocate(Ac, Ac_d)
    deallocate(B, B_d)
    call MPI_FINALIZE(ierr)

contains

    subroutine parse_positive_arg(index, default_value, value, name, rank)
        integer, intent(in) :: index, default_value, rank
        integer, intent(out) :: value
        character(len=*), intent(in) :: name
        character(len=64) :: arg
        integer :: stat, ierr_local

        value = default_value
        if (command_argument_count() >= index) then
            call get_command_argument(index, arg)
            read(arg, *, iostat=stat) value
            if (stat /= 0 .or. value <= 0) then
                if (rank == 0) write(*,*) trim(name), " must be a positive integer."
                call MPI_ABORT(MPI_COMM_WORLD, 1, ierr_local)
            endif
        endif
    end subroutine parse_positive_arg

    subroutine fill_timing_fields(timing_arg, values)
        type(ptdma_timing_cuda), intent(in) :: timing_arg
        real(8), intent(out) :: values(n_timing_fields)

        values(1) = timing_arg%total
        values(2) = timing_arg%local_compute
        values(3) = timing_arg%pack_forward
        values(4) = timing_arg%mpi_forward
        values(5) = timing_arg%unpack_forward
        values(6) = timing_arg%reduced_compute
        values(7) = timing_arg%pack_backward
        values(8) = timing_arg%mpi_backward
        values(9) = timing_arg%unpack_backward
        values(10) = timing_arg%update_compute
        values(11) = timing_arg%local_compute + timing_arg%reduced_compute + timing_arg%update_compute
        values(12) = timing_arg%mpi_forward + timing_arg%mpi_backward
        values(13) = timing_arg%pack_forward + timing_arg%unpack_forward + &
                     timing_arg%pack_backward + timing_arg%unpack_backward
    end subroutine fill_timing_fields

    subroutine write_csv_row(nprocs_arg, n1_arg, n2_arg, n3_arg, nsys_arg, &
                             nrow_min_arg, nrow_max_arg, iter_arg, iterations_arg, &
                             phase_max_arg, total_avg_arg)
        integer, intent(in) :: nprocs_arg, n1_arg, n2_arg, n3_arg
        integer, intent(in) :: nsys_arg, nrow_min_arg, nrow_max_arg
        integer, intent(in) :: iter_arg, iterations_arg
        real(8), intent(in) :: phase_max_arg(n_timing_fields), total_avg_arg
        integer :: i

        write(*,'(A)',advance='no') "tdma,fortran-original,"
        write(*,'(I0,A)',advance='no') nprocs_arg, ","
        write(*,'(I0,A)',advance='no') n1_arg, ","
        write(*,'(I0,A)',advance='no') n2_arg, ","
        write(*,'(I0,A)',advance='no') n3_arg, ","
        write(*,'(I0,A)',advance='no') nsys_arg, ","
        write(*,'(I0,A)',advance='no') nrow_min_arg, ","
        write(*,'(I0,A)',advance='no') nrow_max_arg, ","
        write(*,'(I0,A)',advance='no') iter_arg, ","
        write(*,'(I0,A)',advance='no') iterations_arg, ","
        write(*,'(A)',advance='no') "device,"
        write(*,'(ES24.16,A)',advance='no') phase_max_arg(1), ","
        write(*,'(ES24.16,A)',advance='no') total_avg_arg, ","
        do i = 2, n_timing_fields - 1
            write(*,'(ES24.16,A)',advance='no') phase_max_arg(i), ","
        enddo
        write(*,'(ES24.16)') phase_max_arg(n_timing_fields)
    end subroutine write_csv_row

    subroutine write_correctness_csv_if_requested(solution, n1_arg, n2_arg, n3_arg, &
                                                  nsys_arg, nrow_arg, nrow_min_arg, &
                                                  nrow_max_arg, ia_arg, ib_arg, &
                                                  nprocs_arg, rank_arg)
        integer, intent(in) :: n1_arg, n2_arg, n3_arg, nsys_arg, nrow_arg
        integer, intent(in) :: nrow_min_arg, nrow_max_arg, ia_arg, ib_arg
        integer, intent(in) :: nprocs_arg, rank_arg
        real(8), intent(in) :: solution(0:n1_arg-1,0:n2_arg-1,0:nrow_arg-1)
        real(8), parameter :: expected_value = 1.0d0
        character(len=512) :: output_path
        integer :: path_len, env_status, ierr_local
        integer :: sample_z(3), idx, unit, io_status, file_size
        logical :: file_exists, need_header
        real(8) :: local_sum, local_sumsq, local_linf, local_error
        real(8) :: global_sum, global_sumsq, global_linf, global_error
        real(8) :: local_samples(3), global_samples(3)

        call get_environment_variable("PASCAL_TDMA_CORRECTNESS_OUT", output_path, &
                                      length=path_len, status=env_status)

        local_sum = sum(solution)
        local_sumsq = sum(solution * solution)
        local_linf = maxval(abs(solution))
        local_error = maxval(abs(solution - expected_value))

        sample_z = (/ 0, n3_arg / 2, n3_arg - 1 /)
        local_samples = 0.0d0
        do idx = 1, 3
            if (sample_z(idx) >= ia_arg .and. sample_z(idx) <= ib_arg) then
                local_samples(idx) = solution(0, 0, sample_z(idx) - ia_arg)
            endif
        enddo

        call MPI_REDUCE(local_sum, global_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, ierr_local)
        call MPI_REDUCE(local_sumsq, global_sumsq, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                        0, MPI_COMM_WORLD, ierr_local)
        call MPI_REDUCE(local_linf, global_linf, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                        0, MPI_COMM_WORLD, ierr_local)
        call MPI_REDUCE(local_error, global_error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                        0, MPI_COMM_WORLD, ierr_local)
        call MPI_REDUCE(local_samples, global_samples, 3, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, 0, MPI_COMM_WORLD, ierr_local)

        if (rank_arg /= 0) return
        if (env_status /= 0 .or. path_len <= 0) return

        file_size = 0
        inquire(file=trim(output_path(1:path_len)), exist=file_exists, size=file_size)
        need_header = (.not. file_exists) .or. (file_size == 0)
        open(newunit=unit, file=trim(output_path(1:path_len)), status="unknown", &
             position="append", action="write", iostat=io_status)
        if (io_status /= 0) then
            write(*,*) "Failed to open correctness CSV: ", trim(output_path(1:path_len))
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr_local)
        endif

        if (need_header) call write_correctness_header(unit)
        call write_correctness_row(unit, nprocs_arg, n1_arg, n2_arg, n3_arg, nsys_arg, &
                                   nrow_min_arg, nrow_max_arg, global_sum, global_sumsq, &
                                   global_linf, global_samples, expected_value, global_error)
        close(unit)
    end subroutine write_correctness_csv_if_requested

    subroutine write_correctness_header(unit)
        integer, intent(in) :: unit

        write(unit,'(A)',advance='no') "solver,implementation,nranks,n1,n2,n3,nsys,"
        write(unit,'(A)',advance='no') "nrow_min,nrow_max,mpi_mode,solution_sum,"
        write(unit,'(A)',advance='no') "solution_l2,solution_linf,sample_z0,"
        write(unit,'(A)',advance='no') "sample_zmid,sample_zlast,expected_value,"
        write(unit,'(A)') "max_abs_error_to_expected"
    end subroutine write_correctness_header

    subroutine write_correctness_row(unit, nprocs_arg, n1_arg, n2_arg, n3_arg, &
                                     nsys_arg, nrow_min_arg, nrow_max_arg, &
                                     global_sum_arg, global_sumsq_arg, &
                                     global_linf_arg, global_samples_arg, &
                                     expected_value_arg, global_error_arg)
        integer, intent(in) :: unit, nprocs_arg, n1_arg, n2_arg, n3_arg
        integer, intent(in) :: nsys_arg, nrow_min_arg, nrow_max_arg
        real(8), intent(in) :: global_sum_arg, global_sumsq_arg, global_linf_arg
        real(8), intent(in) :: global_samples_arg(3), expected_value_arg, global_error_arg

        write(unit,'(A)',advance='no') "tdma,fortran-original,"
        write(unit,'(I0,A)',advance='no') nprocs_arg, ","
        write(unit,'(I0,A)',advance='no') n1_arg, ","
        write(unit,'(I0,A)',advance='no') n2_arg, ","
        write(unit,'(I0,A)',advance='no') n3_arg, ","
        write(unit,'(I0,A)',advance='no') nsys_arg, ","
        write(unit,'(I0,A)',advance='no') nrow_min_arg, ","
        write(unit,'(I0,A)',advance='no') nrow_max_arg, ","
        write(unit,'(A)',advance='no') "device,"
        write(unit,'(ES24.16,A)',advance='no') global_sum_arg, ","
        write(unit,'(ES24.16,A)',advance='no') sqrt(global_sumsq_arg), ","
        write(unit,'(ES24.16,A)',advance='no') global_linf_arg, ","
        write(unit,'(ES24.16,A)',advance='no') global_samples_arg(1), ","
        write(unit,'(ES24.16,A)',advance='no') global_samples_arg(2), ","
        write(unit,'(ES24.16,A)',advance='no') global_samples_arg(3), ","
        write(unit,'(ES24.16,A)',advance='no') expected_value_arg, ","
        write(unit,'(ES24.16)') global_error_arg
    end subroutine write_correctness_row

end program example_fortran_profile
