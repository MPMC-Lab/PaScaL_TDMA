#include "pascal_tdma_cuda.hpp"

#include <cuda_runtime.h>
#include <mpi.h>

#include <algorithm>
#include <exception>
#include <iostream>
#include <vector>

namespace {

void ordered_print(const std::vector<double>& d_host,
                   const int nsys,
                   const int nrow,
                   const int n1,
                   const int n2,
                   const int z_first,
                   const int rank,
                   const int nprocs) {
    for (int r = 0; r < nprocs; ++r) {
        if (rank == r) {
            std::cout << "Rank " << rank << " --\n";
            const int first_count = std::min(3, nrow);
            for (int k = 0; k < first_count; ++k) {
                const int sys_a = 0;
                const int sys_b = (n1 / 2) + (n2 / 2) * n1;
                const int sys_c = (n1 - 1) + (n2 - 1) * n1;
                std::cout << "  local k=" << k
                          << " global k=" << (z_first + k)
                          << " : " << d_host[pascal_tdma::index2(sys_a, k, nsys)]
                          << " ... " << d_host[pascal_tdma::index2(sys_b, k, nsys)]
                          << " ... " << d_host[pascal_tdma::index2(sys_c, k, nsys)]
                          << '\n';
            }
            if (nrow > 3) {
                std::cout << "  ...\n";
            }
            for (int k = std::max(first_count, nrow - 2); k < nrow; ++k) {
                const int sys_a = 0;
                const int sys_b = (n1 / 2) + (n2 / 2) * n1;
                const int sys_c = (n1 - 1) + (n2 - 1) * n1;
                std::cout << "  local k=" << k
                          << " global k=" << (z_first + k)
                          << " : " << d_host[pascal_tdma::index2(sys_a, k, nsys)]
                          << " ... " << d_host[pascal_tdma::index2(sys_b, k, nsys)]
                          << " ... " << d_host[pascal_tdma::index2(sys_c, k, nsys)]
                          << '\n';
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

}  // namespace

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank = 0;
    int nprocs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    try {
        int device_count = 0;
        PASCAL_TDMA_CUDA_CHECK(cudaGetDeviceCount(&device_count));
        if (device_count <= 0) {
            throw std::runtime_error("No CUDA device is visible to this MPI rank.");
        }
        const int gpu_rank = rank % device_count;
        PASCAL_TDMA_CUDA_CHECK(cudaSetDevice(gpu_rank));
        PASCAL_TDMA_CUDA_CHECK(cudaDeviceSynchronize());

        const int n1 = 64;
        const int n2 = 64;
        const int n3 = 2048;

        int z_first = 0;
        int z_last = -1;
        pascal_tdma::partition_1d(0, n3 - 1, nprocs, rank, z_first, z_last);
        const int nrow = z_last - z_first + 1;
        const int nsys = n1 * n2;
        const std::size_t nvalues = static_cast<std::size_t>(nsys) * nrow;

        if (rank == 0) {
            std::cout << "MPI ranks: " << nprocs << '\n';
            std::cout << "CUDA devices visible on rank 0: " << device_count << '\n';
            std::cout << "Grid: " << n1 << " x " << n2 << " x " << n3 << '\n';
            std::cout << "MPI mode: "
                      << (pascal_tdma::mpi_mode_from_env() == pascal_tdma::MpiBufferMode::HostStaging
                              ? "host-staging"
                              : "device-direct")
                      << '\n';
        }

        std::vector<double> h_a(nvalues, 1.0);
        std::vector<double> h_b(nvalues, -2.0);
        std::vector<double> h_c(nvalues, 1.0);
        std::vector<double> h_d(nvalues, 0.0);

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

        pascal_tdma::PascalTdmaPlan plan;
        plan.create(nsys, MPI_COMM_WORLD, 128, 128, pascal_tdma::mpi_mode_from_env());
        pascal_tdma::solve(plan, d_a, d_b, d_c, d_d, nsys, nrow);
        PASCAL_TDMA_CUDA_CHECK(cudaDeviceSynchronize());
        plan.destroy();

        PASCAL_TDMA_CUDA_CHECK(cudaMemcpy(h_d.data(), d_d, nvalues * sizeof(double), cudaMemcpyDeviceToHost));
        ordered_print(h_d, nsys, nrow, n1, n2, z_first, rank, nprocs);

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
