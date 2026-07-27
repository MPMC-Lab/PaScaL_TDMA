#include "pascal_tdma_cuda.hpp"

#include <cuda_runtime.h>
#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int kTimingFields = 13;

int parse_positive_int(const int argc,
                       char** argv,
                       const int index,
                       const int default_value,
                       const char* name) {
    if (argc <= index) {
        return default_value;
    }
    char* end = nullptr;
    const long value = std::strtol(argv[index], &end, 10);
    if (end == argv[index] || *end != '\0' || value <= 0) {
        throw std::invalid_argument(std::string(name) + " must be a positive integer");
    }
    return static_cast<int>(value);
}

const char* mpi_mode_name(const pascal_tdma::MpiBufferMode mode) {
    return mode == pascal_tdma::MpiBufferMode::HostStaging ? "host" : "device";
}

std::array<double, kTimingFields> timing_fields(const pascal_tdma::SolveTimings& timings) {
    return {timings.total,
            timings.local_compute,
            timings.pack_forward,
            timings.mpi_forward,
            timings.unpack_forward,
            timings.reduced_compute,
            timings.pack_backward,
            timings.mpi_backward,
            timings.unpack_backward,
            timings.update_compute,
            timings.computation(),
            timings.communication(),
            timings.packing()};
}

void initialize_coefficients(std::vector<double>& h_a,
                             std::vector<double>& h_b,
                             std::vector<double>& h_c,
                             std::vector<double>& h_d,
                             const int nsys,
                             const int nrow,
                             const int rank,
                             const int nprocs) {
    std::fill(h_a.begin(), h_a.end(), 1.0);
    std::fill(h_b.begin(), h_b.end(), -2.0);
    std::fill(h_c.begin(), h_c.end(), 1.0);
    std::fill(h_d.begin(), h_d.end(), 0.0);

    if (rank == 0) {
        for (int sys = 0; sys < nsys; ++sys) {
            h_d[pascal_tdma::index2(sys, 0, nsys)] = -1.0;
        }
    }
    if (rank == nprocs - 1) {
        for (int sys = 0; sys < nsys; ++sys) {
            h_d[pascal_tdma::index2(sys, nrow - 1, nsys)] = -1.0;
        }
    }
}

void write_correctness_csv_if_requested(const std::vector<double>& solution,
                                        const int nsys,
                                        const int nrow,
                                        const int z_first,
                                        const int z_last,
                                        const int nprocs,
                                        const int n1,
                                        const int n2,
                                        const int n3,
                                        const int nrow_min,
                                        const int nrow_max,
                                        const int rank,
                                        const char* mpi_mode) {
    const char* output_path = std::getenv("PASCAL_TDMA_CORRECTNESS_OUT");

    constexpr double expected_value = 1.0;
    double local_sum = 0.0;
    double local_sumsq = 0.0;
    double local_linf = 0.0;
    double local_error = 0.0;
    for (const double value : solution) {
        local_sum += value;
        local_sumsq += value * value;
        local_linf = std::max(local_linf, std::abs(value));
        local_error = std::max(local_error, std::abs(value - expected_value));
    }

    const std::array<int, 3> sample_z = {0, n3 / 2, n3 - 1};
    std::array<double, 3> local_samples = {0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i) {
        if (sample_z[i] >= z_first && sample_z[i] <= z_last) {
            const int local_row = sample_z[i] - z_first;
            local_samples[i] = solution[pascal_tdma::index2(0, local_row, nsys)];
        }
    }

    double global_sum = 0.0;
    double global_sumsq = 0.0;
    double global_linf = 0.0;
    double global_error = 0.0;
    std::array<double, 3> global_samples = {0.0, 0.0, 0.0};
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_sumsq, &global_sumsq, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_linf, &global_linf, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(local_samples.data(), global_samples.data(), 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        return;
    }
    if (output_path == nullptr || output_path[0] == '\0') {
        return;
    }

    std::ifstream existing(output_path);
    const bool write_header = !existing.good() ||
                              existing.peek() == std::ifstream::traits_type::eof();
    existing.close();

    std::ofstream out(output_path, std::ios::app);
    if (!out) {
        throw std::runtime_error(std::string("Failed to open correctness CSV: ") + output_path);
    }

    if (write_header) {
        out << "solver,implementation,nranks,n1,n2,n3,nsys,nrow_min,nrow_max,"
            << "mpi_mode,solution_sum,solution_l2,solution_linf,sample_z0,"
            << "sample_zmid,sample_zlast,expected_value,max_abs_error_to_expected\n";
    }

    out << std::setprecision(16)
        << "tdma,cuda-cxx"
        << ',' << nprocs
        << ',' << n1
        << ',' << n2
        << ',' << n3
        << ',' << nsys
        << ',' << nrow_min
        << ',' << nrow_max
        << ',' << mpi_mode
        << ',' << global_sum
        << ',' << std::sqrt(global_sumsq)
        << ',' << global_linf
        << ',' << global_samples[0]
        << ',' << global_samples[1]
        << ',' << global_samples[2]
        << ',' << expected_value
        << ',' << global_error
        << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank = 0;
    int nprocs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    try {
        const int n1 = parse_positive_int(argc, argv, 1, 64, "n1");
        const int n2 = parse_positive_int(argc, argv, 2, 64, "n2");
        const int n3 = parse_positive_int(argc, argv, 3, 2048, "n3");
        const int iterations = parse_positive_int(argc, argv, 4, 10, "iterations");
        const int tdma_threads = parse_positive_int(argc, argv, 5, 128, "tdma_threads");
        const int reduced_threads = parse_positive_int(argc, argv, 6, 128, "reduced_threads");
        if (argc > 7) {
            throw std::invalid_argument(
                "usage: example_cuda_cxx_profile [n1] [n2] [n3] [iterations] [tdma_threads] [reduced_threads]");
        }

        int device_count = 0;
        PASCAL_TDMA_CUDA_CHECK(cudaGetDeviceCount(&device_count));
        if (device_count <= 0) {
            throw std::runtime_error("No CUDA device is visible to this MPI rank.");
        }
        PASCAL_TDMA_CUDA_CHECK(cudaSetDevice(rank % device_count));
        PASCAL_TDMA_CUDA_CHECK(cudaDeviceSynchronize());

        int z_first = 0;
        int z_last = -1;
        pascal_tdma::partition_1d(0, n3 - 1, nprocs, rank, z_first, z_last);
        const int nrow = z_last - z_first + 1;
        const int nsys = n1 * n2;
        const std::size_t nvalues = static_cast<std::size_t>(nsys) * nrow;

        int nrow_min = 0;
        int nrow_max = 0;
        MPI_Reduce(&nrow, &nrow_min, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&nrow, &nrow_max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

        std::vector<double> h_a(nvalues);
        std::vector<double> h_b(nvalues);
        std::vector<double> h_c(nvalues);
        std::vector<double> h_d(nvalues);
        initialize_coefficients(h_a, h_b, h_c, h_d, nsys, nrow, rank, nprocs);

        double* d_a = nullptr;
        double* d_b = nullptr;
        double* d_c = nullptr;
        double* d_d = nullptr;
        PASCAL_TDMA_CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_a), nvalues * sizeof(double)));
        PASCAL_TDMA_CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_b), nvalues * sizeof(double)));
        PASCAL_TDMA_CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_c), nvalues * sizeof(double)));
        PASCAL_TDMA_CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&d_d), nvalues * sizeof(double)));

        PASCAL_TDMA_CUDA_CHECK(cudaMemcpy(d_a, h_a.data(), nvalues * sizeof(double), cudaMemcpyHostToDevice));
        PASCAL_TDMA_CUDA_CHECK(cudaMemcpy(d_b, h_b.data(), nvalues * sizeof(double), cudaMemcpyHostToDevice));
        PASCAL_TDMA_CUDA_CHECK(cudaMemcpy(d_c, h_c.data(), nvalues * sizeof(double), cudaMemcpyHostToDevice));
        PASCAL_TDMA_CUDA_CHECK(cudaMemcpy(d_d, h_d.data(), nvalues * sizeof(double), cudaMemcpyHostToDevice));

        const pascal_tdma::MpiBufferMode mpi_mode = pascal_tdma::mpi_mode_from_env();
        pascal_tdma::PascalTdmaPlan plan;
        plan.create(nsys, MPI_COMM_WORLD, tdma_threads, reduced_threads, mpi_mode);

        if (rank == 0) {
            std::cout << "solver,implementation,nranks,n1,n2,n3,nsys,nrow_min,nrow_max,"
                      << "iter,iterations,mpi_mode,total_s_max,total_s_avg,"
                      << "local_compute_s_max,pack_forward_s_max,mpi_forward_s_max,"
                      << "unpack_forward_s_max,reduced_compute_s_max,pack_backward_s_max,"
                      << "mpi_backward_s_max,unpack_backward_s_max,update_compute_s_max,"
                      << "compute_s_max,communication_s_max,packing_s_max\n";
        }

        for (int iter = 0; iter < iterations; ++iter) {
            MPI_Barrier(MPI_COMM_WORLD);

            pascal_tdma::SolveTimings timings;
            pascal_tdma::solve_profiled(plan, d_a, d_b, d_c, d_d, nsys, nrow, &timings);

            if (iter == 0) {
                PASCAL_TDMA_CUDA_CHECK(cudaMemcpy(h_d.data(), d_d, nvalues * sizeof(double),
                                                  cudaMemcpyDeviceToHost));
                write_correctness_csv_if_requested(h_d, nsys, nrow, z_first, z_last, nprocs,
                                                   n1, n2, n3, nrow_min, nrow_max, rank,
                                                   mpi_mode_name(mpi_mode));
            }

            const std::array<double, kTimingFields> local_fields = timing_fields(timings);
            std::array<double, kTimingFields> max_fields{};
            std::array<double, kTimingFields> sum_fields{};
            MPI_Reduce(local_fields.data(), max_fields.data(), kTimingFields,
                       MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(local_fields.data(), sum_fields.data(), kTimingFields,
                       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            if (rank == 0) {
                std::cout << std::setprecision(12)
                          << "tdma,cuda-cxx"
                          << ',' << nprocs
                          << ',' << n1
                          << ',' << n2
                          << ',' << n3
                          << ',' << nsys
                          << ',' << nrow_min
                          << ',' << nrow_max
                          << ',' << iter
                          << ',' << iterations
                          << ',' << mpi_mode_name(mpi_mode)
                          << ',' << max_fields[0]
                          << ',' << (sum_fields[0] / nprocs);
                for (int field = 1; field < kTimingFields; ++field) {
                    std::cout << ',' << max_fields[field];
                }
                std::cout << '\n';
            }
        }

        plan.destroy();
        cudaFree(d_a);
        cudaFree(d_b);
        cudaFree(d_c);
        cudaFree(d_d);
    } catch (const std::exception& e) {
        std::cerr << "Rank " << rank << " error: " << e.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Finalize();
    return 0;
}
