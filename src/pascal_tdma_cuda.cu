#include "pascal_tdma_cuda.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <sstream>
#include <utility>

namespace pascal_tdma {
namespace {

inline int ceil_div(const int a, const int b) {
    return (a + b - 1) / b;
}

inline int desc_at(const std::vector<int>& values, const int dim, const int rank) {
    return values[2 * rank + dim];
}

inline std::size_t total_count(const std::vector<int>& values) {
    return static_cast<std::size_t>(
        std::accumulate(values.begin(), values.end(), 0));
}

__device__ __forceinline__ std::size_t d_index2(const int sys, const int row, const int nsys) {
    return static_cast<std::size_t>(sys) + static_cast<std::size_t>(row) * nsys;
}

__global__ void tdma_many_kernel(double* a,
                                 double* b,
                                 double* c,
                                 double* d,
                                 const int nsys,
                                 const int nrow) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= nsys || nrow <= 0) {
        return;
    }

    double b1 = b[d_index2(i, 0, nsys)];
    double c1 = c[d_index2(i, 0, nsys)];
    double d1 = d[d_index2(i, 0, nsys)];

    d1 /= b1;
    c1 /= b1;

    d[d_index2(i, 0, nsys)] = d1;
    c[d_index2(i, 0, nsys)] = c1;

    for (int j = 1; j < nrow; ++j) {
        const double c0 = c1;
        const double d0 = d1;

        const double a1 = a[d_index2(i, j, nsys)];
        b1 = b[d_index2(i, j, nsys)];
        c1 = c[d_index2(i, j, nsys)];
        d1 = d[d_index2(i, j, nsys)];

        const double r = 1.0 / (b1 - a1 * c0);
        d1 = r * (d1 - a1 * d0);
        c1 = r * c1;

        d[d_index2(i, j, nsys)] = d1;
        c[d_index2(i, j, nsys)] = c1;
    }

    for (int j = nrow - 2; j >= 0; --j) {
        const double c0 = c[d_index2(i, j, nsys)];
        double d0 = d[d_index2(i, j, nsys)];
        d0 -= c0 * d1;
        d1 = d0;
        d[d_index2(i, j, nsys)] = d0;
    }
}

__global__ void tdma_modified_kernel(double* a,
                                     double* b,
                                     double* c,
                                     double* d,
                                     double* a_rd,
                                     double* b_rd,
                                     double* c_rd,
                                     double* d_rd,
                                     const int nsys,
                                     const int nrow) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= nsys) {
        return;
    }

    double a0 = a[d_index2(i, 0, nsys)];
    double b0 = b[d_index2(i, 0, nsys)];
    double c0 = c[d_index2(i, 0, nsys)];
    double d0 = d[d_index2(i, 0, nsys)];

    a0 /= b0;
    c0 /= b0;
    d0 /= b0;

    a[d_index2(i, 0, nsys)] = a0;
    c[d_index2(i, 0, nsys)] = c0;
    d[d_index2(i, 0, nsys)] = d0;

    double a1 = a[d_index2(i, 1, nsys)];
    double b1 = b[d_index2(i, 1, nsys)];
    double c1 = c[d_index2(i, 1, nsys)];
    double d1 = d[d_index2(i, 1, nsys)];

    a1 /= b1;
    c1 /= b1;
    d1 /= b1;

    a[d_index2(i, 1, nsys)] = a1;
    c[d_index2(i, 1, nsys)] = c1;
    d[d_index2(i, 1, nsys)] = d1;

    for (int j = 2; j < nrow; ++j) {
        a0 = a1;
        c0 = c1;
        d0 = d1;

        a1 = a[d_index2(i, j, nsys)];
        b1 = b[d_index2(i, j, nsys)];
        c1 = c[d_index2(i, j, nsys)];
        d1 = d[d_index2(i, j, nsys)];

        const double r = 1.0 / (b1 - a1 * c0);
        d1 = r * (d1 - a1 * d0);
        c1 = r * c1;
        a1 = -r * a1 * a0;

        a[d_index2(i, j, nsys)] = a1;
        c[d_index2(i, j, nsys)] = c1;
        d[d_index2(i, j, nsys)] = d1;
    }

    a_rd[d_index2(i, 1, nsys)] = a1;
    b_rd[d_index2(i, 1, nsys)] = 1.0;
    c_rd[d_index2(i, 1, nsys)] = c1;
    d_rd[d_index2(i, 1, nsys)] = d1;

    a1 = a0;
    c1 = c0;
    d1 = d0;

    for (int j = nrow - 3; j >= 1; --j) {
        a0 = a[d_index2(i, j, nsys)];
        d0 = d[d_index2(i, j, nsys)];
        c0 = c[d_index2(i, j, nsys)];

        a0 = a0 - c0 * a1;
        d0 = d0 - c0 * d1;
        c0 = -c0 * c1;

        a1 = a0;
        d1 = d0;
        c1 = c0;

        a[d_index2(i, j, nsys)] = a0;
        d[d_index2(i, j, nsys)] = d0;
        c[d_index2(i, j, nsys)] = c0;
    }

    a0 = a[d_index2(i, 0, nsys)];
    d0 = d[d_index2(i, 0, nsys)];
    c0 = c[d_index2(i, 0, nsys)];

    const double r = 1.0 / (1.0 - a1 * c0);
    a0 = r * a0;
    d0 = r * (d0 - c0 * d1);
    c0 = -r * c0 * c1;

    a[d_index2(i, 0, nsys)] = a0;
    d[d_index2(i, 0, nsys)] = d0;
    c[d_index2(i, 0, nsys)] = c0;

    a_rd[d_index2(i, 0, nsys)] = a0;
    b_rd[d_index2(i, 0, nsys)] = 1.0;
    c_rd[d_index2(i, 0, nsys)] = c0;
    d_rd[d_index2(i, 0, nsys)] = d0;
}

__global__ void pascal_update_kernel(double* a,
                                     double* b,
                                     double* c,
                                     double* d,
                                     double* d_rd,
                                     const int nsys,
                                     const int nrow) {
    (void)b;
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= nsys) {
        return;
    }

    const double ds = d_rd[d_index2(i, 0, nsys)];
    const double de = d_rd[d_index2(i, 1, nsys)];

    d[d_index2(i, 0, nsys)] = ds;
    d[d_index2(i, nrow - 1, nsys)] = de;

    for (int j = 1; j < nrow - 1; ++j) {
        d[d_index2(i, j, nsys)] -=
            a[d_index2(i, j, nsys)] * ds + c[d_index2(i, j, nsys)] * de;
    }
}

__global__ void pack_2d_kernel(const double* a,
                               const int n1,
                               const int n2,
                               const int sub0,
                               const int sub1,
                               const int start0,
                               const int start1,
                               double* buffer,
                               const std::size_t buffer_size,
                               const int buffer_offset) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    const int j = blockIdx.y;
    const int index_i = i + start0;
    const int index_j = j + start1;
    if (i < sub0 && j < sub1 && index_i < n1 && index_j < n2) {
        const std::size_t out = static_cast<std::size_t>(buffer_offset) + i + j * sub0;
        if (out < buffer_size) {
            buffer[out] = a[d_index2(index_i, index_j, n1)];
        }
    }
}

__global__ void unpack_2d_kernel(double* a,
                                 const int n1,
                                 const int n2,
                                 const int sub0,
                                 const int sub1,
                                 const int start0,
                                 const int start1,
                                 const double* buffer,
                                 const std::size_t buffer_size,
                                 const int buffer_offset) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    const int j = blockIdx.y;
    const int index_i = i + start0;
    const int index_j = j + start1;
    if (i < sub0 && j < sub1 && index_i < n1 && index_j < n2) {
        const std::size_t in = static_cast<std::size_t>(buffer_offset) + i + j * sub0;
        if (in < buffer_size) {
            a[d_index2(index_i, index_j, n1)] = buffer[in];
        }
    }
}

__global__ void init_2d_kernel(double* a, const int n1, const int n2) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    const int j = blockIdx.y;
    if (i < n1 && j < n2) {
        a[d_index2(i, j, n1)] = 1.0;
    }
}

void launch_tdma_many(double* a,
                      double* b,
                      double* c,
                      double* d,
                      const int nsys,
                      const int nrow,
                      const dim3 blocks,
                      const dim3 threads,
                      cudaStream_t stream) {
    tdma_many_kernel<<<blocks, threads, 0, stream>>>(a, b, c, d, nsys, nrow);
    PASCAL_TDMA_CUDA_CHECK(cudaGetLastError());
}

void launch_pack(const double* src,
                 const int n1,
                 const int n2,
                 const int sub0,
                 const int sub1,
                 const int start0,
                 const int start1,
                 double* buffer,
                 const std::size_t buffer_size,
                 const int buffer_offset,
                 const dim3 threads,
                 cudaStream_t stream) {
    if (sub0 == 0 || sub1 == 0) {
        return;
    }
    const dim3 blocks(ceil_div(sub0, static_cast<int>(threads.x)), sub1, 1);
    pack_2d_kernel<<<blocks, threads, 0, stream>>>(
        src, n1, n2, sub0, sub1, start0, start1, buffer, buffer_size, buffer_offset);
    PASCAL_TDMA_CUDA_CHECK(cudaGetLastError());
}

void launch_unpack(double* dst,
                   const int n1,
                   const int n2,
                   const int sub0,
                   const int sub1,
                   const int start0,
                   const int start1,
                   const double* buffer,
                   const std::size_t buffer_size,
                   const int buffer_offset,
                   const dim3 threads,
                   cudaStream_t stream) {
    if (sub0 == 0 || sub1 == 0) {
        return;
    }
    const dim3 blocks(ceil_div(sub0, static_cast<int>(threads.x)), sub1, 1);
    unpack_2d_kernel<<<blocks, threads, 0, stream>>>(
        dst, n1, n2, sub0, sub1, start0, start1, buffer, buffer_size, buffer_offset);
    PASCAL_TDMA_CUDA_CHECK(cudaGetLastError());
}

void alltoallv_double(double* send_dev,
                      const std::size_t send_size,
                      const std::vector<int>& send_counts,
                      const std::vector<int>& send_displs,
                      double* recv_dev,
                      const std::size_t recv_size,
                      const std::vector<int>& recv_counts,
                      const std::vector<int>& recv_displs,
                      MPI_Comm comm,
                      const MpiBufferMode mode,
                      cudaStream_t stream) {
    PASCAL_TDMA_CUDA_CHECK(cudaStreamSynchronize(stream));

    int ierr = MPI_SUCCESS;
    if (mode == MpiBufferMode::DeviceDirect) {
        ierr = MPI_Alltoallv(send_dev,
                             send_counts.data(),
                             send_displs.data(),
                             MPI_DOUBLE,
                             recv_dev,
                             recv_counts.data(),
                             recv_displs.data(),
                             MPI_DOUBLE,
                             comm);
    } else {
        std::vector<double> send_host(send_size);
        std::vector<double> recv_host(recv_size);
        PASCAL_TDMA_CUDA_CHECK(cudaMemcpyAsync(send_host.data(),
                                               send_dev,
                                               send_size * sizeof(double),
                                               cudaMemcpyDeviceToHost,
                                               stream));
        PASCAL_TDMA_CUDA_CHECK(cudaStreamSynchronize(stream));
        ierr = MPI_Alltoallv(send_host.data(),
                             send_counts.data(),
                             send_displs.data(),
                             MPI_DOUBLE,
                             recv_host.data(),
                             recv_counts.data(),
                             recv_displs.data(),
                             MPI_DOUBLE,
                             comm);
        if (ierr == MPI_SUCCESS) {
            PASCAL_TDMA_CUDA_CHECK(cudaMemcpyAsync(recv_dev,
                                                   recv_host.data(),
                                                   recv_size * sizeof(double),
                                                   cudaMemcpyHostToDevice,
                                                   stream));
        }
    }

    if (ierr != MPI_SUCCESS) {
        throw std::runtime_error("MPI_Alltoallv failed in pascal_tdma::alltoallv_double");
    }
}

void swap_plan(PascalTdmaPlan& lhs, PascalTdmaPlan& rhs) noexcept {
    using std::swap;
    swap(lhs.comm, rhs.comm);
    swap(lhs.rank, rhs.rank);
    swap(lhs.nprocs, rhs.nprocs);
    swap(lhs.nsys, rhs.nsys);
    swap(lhs.local_nsys, rhs.local_nsys);
    swap(lhs.local_first_sys, rhs.local_first_sys);
    swap(lhs.local_last_sys, rhs.local_last_sys);
    swap(lhs.nrd_global, rhs.nrd_global);
    swap(lhs.nrd_local, rhs.nrd_local);
    swap(lhs.ntr_global, rhs.ntr_global);
    swap(lhs.ntr_local, rhs.ntr_local);
    swap(lhs.gather_nrd_local, rhs.gather_nrd_local);
    swap(lhs.gather_ntr_local, rhs.gather_ntr_local);
    swap(lhs.gather_nrd_start, rhs.gather_nrd_start);
    swap(lhs.gather_ntr_start, rhs.gather_ntr_start);
    swap(lhs.counts_a, rhs.counts_a);
    swap(lhs.displs_a, rhs.displs_a);
    swap(lhs.counts_b, rhs.counts_b);
    swap(lhs.displs_b, rhs.displs_b);
    swap(lhs.big_counts_a, rhs.big_counts_a);
    swap(lhs.big_displs_a, rhs.big_displs_a);
    swap(lhs.big_counts_b, rhs.big_counts_b);
    swap(lhs.big_displs_b, rhs.big_displs_b);
    swap(lhs.ard, rhs.ard);
    swap(lhs.brd, rhs.brd);
    swap(lhs.crd, rhs.crd);
    swap(lhs.drd, rhs.drd);
    swap(lhs.atr, rhs.atr);
    swap(lhs.btr, rhs.btr);
    swap(lhs.ctr, rhs.ctr);
    swap(lhs.dtr, rhs.dtr);
    swap(lhs.bigbuf_a, rhs.bigbuf_a);
    swap(lhs.bigbuf_b, rhs.bigbuf_b);
    swap(lhs.bigbuf_a_size, rhs.bigbuf_a_size);
    swap(lhs.bigbuf_b_size, rhs.bigbuf_b_size);
    swap(lhs.buf_a_size, rhs.buf_a_size);
    swap(lhs.buf_b_size, rhs.buf_b_size);
    swap(lhs.threads_tdma, rhs.threads_tdma);
    swap(lhs.blocks_tdma, rhs.blocks_tdma);
    swap(lhs.threads_reduced, rhs.threads_reduced);
    swap(lhs.blocks_reduced, rhs.blocks_reduced);
    swap(lhs.threads_pack, rhs.threads_pack);
    swap(lhs.threads_init, rhs.threads_init);
    swap(lhs.mpi_mode, rhs.mpi_mode);
    swap(lhs.created, rhs.created);
}

}  // namespace

void cuda_check(const cudaError_t status, const char* expr, const char* file, const int line) {
    if (status != cudaSuccess) {
        std::ostringstream oss;
        oss << file << ':' << line << ": CUDA call failed: " << expr
            << " -> " << cudaGetErrorString(status);
        throw CudaError(oss.str());
    }
}

void partition_1d(const int start,
                  const int end,
                  const int nprocs,
                  const int rank,
                  int& first,
                  int& last) {
    const int n = end - start + 1;
    const int base = n / nprocs;
    const int rem = n % nprocs;
    first = rank * base + start + std::min(rank, rem);
    last = first + base - 1;
    if (rem > rank) {
        ++last;
    }
}

MpiBufferMode mpi_mode_from_env() {
    const char* value = std::getenv("PASCAL_TDMA_MPI_MODE");
    if (value != nullptr && std::strcmp(value, "host") == 0) {
        return MpiBufferMode::HostStaging;
    }
    return MpiBufferMode::DeviceDirect;
}

PascalTdmaPlan::PascalTdmaPlan(PascalTdmaPlan&& other) noexcept {
    swap_plan(*this, other);
}

PascalTdmaPlan& PascalTdmaPlan::operator=(PascalTdmaPlan&& other) noexcept {
    if (this != &other) {
        destroy();
        swap_plan(*this, other);
    }
    return *this;
}

PascalTdmaPlan::~PascalTdmaPlan() {
    destroy();
}

void PascalTdmaPlan::create(const int nsys_in,
                            MPI_Comm comm_in,
                            int tdma_threads,
                            int reduced_threads,
                            const MpiBufferMode mode) {
    destroy();
    if (nsys_in <= 0) {
        throw std::invalid_argument("PascalTdmaPlan::create requires nsys > 0");
    }
    if (tdma_threads <= 0) {
        tdma_threads = 128;
    }
    if (reduced_threads <= 0) {
        reduced_threads = 128;
    }

    MPI_Comm_dup(comm_in, &comm);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);
    nsys = nsys_in;
    mpi_mode = mode;

    partition_1d(0, nsys - 1, nprocs, rank, local_first_sys, local_last_sys);
    local_nsys = local_last_sys - local_first_sys + 1;

    nrd_global[0] = nsys;
    nrd_global[1] = 2;
    nrd_local[0] = local_nsys;
    nrd_local[1] = 2;
    ntr_global[0] = local_nsys;
    ntr_global[1] = 2 * nprocs;
    ntr_local[0] = local_nsys;
    ntr_local[1] = 2;

    gather_nrd_local.assign(2 * nprocs, 0);
    gather_ntr_local.assign(2 * nprocs, 0);
    gather_nrd_start.assign(2 * nprocs, 0);
    gather_ntr_start.assign(2 * nprocs, 0);

    counts_a.assign(nprocs, 0);
    displs_a.assign(nprocs, 0);
    counts_b.assign(nprocs, 0);
    displs_b.assign(nprocs, 0);
    big_counts_a.assign(nprocs, 0);
    big_displs_a.assign(nprocs, 0);
    big_counts_b.assign(nprocs, 0);
    big_displs_b.assign(nprocs, 0);

    int prefix_rows = 0;
    int prefix_a = 0;
    int prefix_b = 0;
    for (int r = 0; r < nprocs; ++r) {
        int first = 0;
        int last = -1;
        partition_1d(0, nsys - 1, nprocs, r, first, last);
        const int rows = last - first + 1;

        gather_nrd_local[2 * r + 0] = rows;
        gather_nrd_local[2 * r + 1] = 2;
        gather_ntr_local[2 * r + 0] = local_nsys;
        gather_ntr_local[2 * r + 1] = 2;

        gather_nrd_start[2 * r + 0] = prefix_rows;
        gather_nrd_start[2 * r + 1] = 0;
        gather_ntr_start[2 * r + 0] = 0;
        gather_ntr_start[2 * r + 1] = 2 * r;

        counts_a[r] = rows * 2;
        displs_a[r] = prefix_a;
        prefix_a += counts_a[r];

        counts_b[r] = local_nsys * 2;
        displs_b[r] = prefix_b;
        prefix_b += counts_b[r];

        prefix_rows += rows;
    }

    for (int r = 0; r < nprocs; ++r) {
        big_counts_a[r] = 3 * counts_a[r];
        big_displs_a[r] = 3 * displs_a[r];
        big_counts_b[r] = 3 * counts_b[r];
        big_displs_b[r] = 3 * displs_b[r];
    }

    buf_a_size = total_count(counts_a);
    buf_b_size = total_count(counts_b);
    bigbuf_a_size = total_count(big_counts_a);
    bigbuf_b_size = total_count(big_counts_b);

    PASCAL_TDMA_CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&ard), sizeof(double) * nsys * 2));
    PASCAL_TDMA_CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&brd), sizeof(double) * nsys * 2));
    PASCAL_TDMA_CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&crd), sizeof(double) * nsys * 2));
    PASCAL_TDMA_CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&drd), sizeof(double) * nsys * 2));

    const std::size_t tr_size = static_cast<std::size_t>(local_nsys) * 2 * nprocs;
    PASCAL_TDMA_CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&atr), sizeof(double) * tr_size));
    PASCAL_TDMA_CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&btr), sizeof(double) * tr_size));
    PASCAL_TDMA_CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&ctr), sizeof(double) * tr_size));
    PASCAL_TDMA_CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&dtr), sizeof(double) * tr_size));

    PASCAL_TDMA_CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&bigbuf_a), sizeof(double) * bigbuf_a_size));
    PASCAL_TDMA_CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&bigbuf_b), sizeof(double) * bigbuf_b_size));

    threads_tdma = dim3(static_cast<unsigned int>(tdma_threads), 1, 1);
    blocks_tdma = dim3(static_cast<unsigned int>(ceil_div(nsys, tdma_threads)), 1, 1);
    threads_reduced = dim3(static_cast<unsigned int>(reduced_threads), 1, 1);
    blocks_reduced = dim3(static_cast<unsigned int>(ceil_div(local_nsys, reduced_threads)), 1, 1);

    const int max_rows = *std::max_element(counts_a.begin(), counts_a.end()) / 2;
    const int pack_threads = std::min(128, std::max(1, max_rows));
    threads_pack = dim3(static_cast<unsigned int>(pack_threads), 1, 1);
    threads_init = dim3(static_cast<unsigned int>(reduced_threads), 1, 1);

    created = true;
}

void PascalTdmaPlan::destroy() noexcept {
    if (ard != nullptr) {
        cudaFree(ard);
    }
    if (brd != nullptr) {
        cudaFree(brd);
    }
    if (crd != nullptr) {
        cudaFree(crd);
    }
    if (drd != nullptr) {
        cudaFree(drd);
    }
    if (atr != nullptr) {
        cudaFree(atr);
    }
    if (btr != nullptr) {
        cudaFree(btr);
    }
    if (ctr != nullptr) {
        cudaFree(ctr);
    }
    if (dtr != nullptr) {
        cudaFree(dtr);
    }
    if (bigbuf_a != nullptr) {
        cudaFree(bigbuf_a);
    }
    if (bigbuf_b != nullptr) {
        cudaFree(bigbuf_b);
    }

    ard = brd = crd = drd = nullptr;
    atr = btr = ctr = dtr = nullptr;
    bigbuf_a = bigbuf_b = nullptr;

    int initialized = 0;
    int finalized = 0;
    MPI_Initialized(&initialized);
    if (initialized) {
        MPI_Finalized(&finalized);
    }
    if (initialized && !finalized && comm != MPI_COMM_NULL) {
        MPI_Comm_free(&comm);
    }
    comm = MPI_COMM_NULL;
    created = false;
}

namespace {

double wall_time() {
    return MPI_Wtime();
}

void add_elapsed(SolveTimings* timings,
                 double SolveTimings::*field,
                 const double start_time) {
    if (timings != nullptr) {
        timings->*field += wall_time() - start_time;
    }
}

void sync_for_timing(SolveTimings* timings, cudaStream_t stream) {
    if (timings != nullptr) {
        PASCAL_TDMA_CUDA_CHECK(cudaStreamSynchronize(stream));
    }
}

void solve_impl(PascalTdmaPlan& plan,
                double* a_dev,
                double* b_dev,
                double* c_dev,
                double* d_dev,
                const int nsys,
                const int nrow,
                SolveTimings* timings,
                cudaStream_t stream) {
    const double total_start = wall_time();

    if (!plan.created) {
        throw std::runtime_error("pascal_tdma::solve called with an uncreated plan");
    }
    if (plan.nsys != nsys) {
        throw std::invalid_argument("pascal_tdma::solve nsys does not match the plan");
    }
    if (nrow <= 0) {
        throw std::invalid_argument("pascal_tdma::solve requires nrow > 0");
    }

    if (plan.nprocs == 1) {
        const double local_start = wall_time();
        launch_tdma_many(a_dev, b_dev, c_dev, d_dev, nsys, nrow,
                         plan.blocks_tdma, plan.threads_tdma, stream);
        sync_for_timing(timings, stream);
        add_elapsed(timings, &SolveTimings::local_compute, local_start);
        add_elapsed(timings, &SolveTimings::total, total_start);
        return;
    }

    if (nrow < 2) {
        throw std::invalid_argument("parallel modified TDMA requires nrow >= 2");
    }

    const double local_start = wall_time();
    tdma_modified_kernel<<<plan.blocks_tdma, plan.threads_tdma, 0, stream>>>(
        a_dev, b_dev, c_dev, d_dev, plan.ard, plan.brd, plan.crd, plan.drd, nsys, nrow);
    PASCAL_TDMA_CUDA_CHECK(cudaGetLastError());
    sync_for_timing(timings, stream);
    add_elapsed(timings, &SolveTimings::local_compute, local_start);

    const double pack_forward_start = wall_time();
    for (int r = 0; r < plan.nprocs; ++r) {
        const int sub0 = desc_at(plan.gather_nrd_local, 0, r);
        const int sub1 = desc_at(plan.gather_nrd_local, 1, r);
        const int start0 = desc_at(plan.gather_nrd_start, 0, r);
        const int start1 = desc_at(plan.gather_nrd_start, 1, r);
        launch_pack(plan.ard, plan.nrd_global[0], plan.nrd_global[1],
                    sub0, sub1, start0, start1, plan.bigbuf_a, plan.bigbuf_a_size,
                    plan.big_displs_a[r] + 0 * plan.counts_a[r], plan.threads_pack, stream);
        launch_pack(plan.crd, plan.nrd_global[0], plan.nrd_global[1],
                    sub0, sub1, start0, start1, plan.bigbuf_a, plan.bigbuf_a_size,
                    plan.big_displs_a[r] + 1 * plan.counts_a[r], plan.threads_pack, stream);
        launch_pack(plan.drd, plan.nrd_global[0], plan.nrd_global[1],
                    sub0, sub1, start0, start1, plan.bigbuf_a, plan.bigbuf_a_size,
                    plan.big_displs_a[r] + 2 * plan.counts_a[r], plan.threads_pack, stream);
    }
    sync_for_timing(timings, stream);
    add_elapsed(timings, &SolveTimings::pack_forward, pack_forward_start);

    const double mpi_forward_start = wall_time();
    alltoallv_double(plan.bigbuf_a, plan.bigbuf_a_size,
                     plan.big_counts_a, plan.big_displs_a,
                     plan.bigbuf_b, plan.bigbuf_b_size,
                     plan.big_counts_b, plan.big_displs_b,
                     plan.comm, plan.mpi_mode, stream);
    sync_for_timing(timings, stream);
    add_elapsed(timings, &SolveTimings::mpi_forward, mpi_forward_start);

    const double unpack_forward_start = wall_time();
    for (int r = 0; r < plan.nprocs; ++r) {
        const int sub0 = desc_at(plan.gather_ntr_local, 0, r);
        const int sub1 = desc_at(plan.gather_ntr_local, 1, r);
        const int start0 = desc_at(plan.gather_ntr_start, 0, r);
        const int start1 = desc_at(plan.gather_ntr_start, 1, r);
        launch_unpack(plan.atr, plan.ntr_global[0], plan.ntr_global[1],
                      sub0, sub1, start0, start1, plan.bigbuf_b, plan.bigbuf_b_size,
                      plan.big_displs_b[r] + 0 * plan.counts_b[r], plan.threads_pack, stream);
        launch_unpack(plan.ctr, plan.ntr_global[0], plan.ntr_global[1],
                      sub0, sub1, start0, start1, plan.bigbuf_b, plan.bigbuf_b_size,
                      plan.big_displs_b[r] + 1 * plan.counts_b[r], plan.threads_pack, stream);
        launch_unpack(plan.dtr, plan.ntr_global[0], plan.ntr_global[1],
                      sub0, sub1, start0, start1, plan.bigbuf_b, plan.bigbuf_b_size,
                      plan.big_displs_b[r] + 2 * plan.counts_b[r], plan.threads_pack, stream);
    }
    sync_for_timing(timings, stream);
    add_elapsed(timings, &SolveTimings::unpack_forward, unpack_forward_start);

    const double reduced_start = wall_time();
    const dim3 init_blocks(ceil_div(plan.ntr_global[0], static_cast<int>(plan.threads_init.x)),
                           plan.ntr_global[1],
                           1);
    init_2d_kernel<<<init_blocks, plan.threads_init, 0, stream>>>(
        plan.btr, plan.ntr_global[0], plan.ntr_global[1]);
    PASCAL_TDMA_CUDA_CHECK(cudaGetLastError());

    launch_tdma_many(plan.atr, plan.btr, plan.ctr, plan.dtr,
                     plan.ntr_global[0], plan.ntr_global[1],
                     plan.blocks_reduced, plan.threads_reduced, stream);
    sync_for_timing(timings, stream);
    add_elapsed(timings, &SolveTimings::reduced_compute, reduced_start);

    const double pack_backward_start = wall_time();
    for (int r = 0; r < plan.nprocs; ++r) {
        const int sub0 = desc_at(plan.gather_ntr_local, 0, r);
        const int sub1 = desc_at(plan.gather_ntr_local, 1, r);
        const int start0 = desc_at(plan.gather_ntr_start, 0, r);
        const int start1 = desc_at(plan.gather_ntr_start, 1, r);
        launch_pack(plan.dtr, plan.ntr_global[0], plan.ntr_global[1],
                    sub0, sub1, start0, start1, plan.bigbuf_b, plan.buf_b_size,
                    plan.displs_b[r], plan.threads_pack, stream);
    }
    sync_for_timing(timings, stream);
    add_elapsed(timings, &SolveTimings::pack_backward, pack_backward_start);

    const double mpi_backward_start = wall_time();
    alltoallv_double(plan.bigbuf_b, plan.buf_b_size,
                     plan.counts_b, plan.displs_b,
                     plan.bigbuf_a, plan.buf_a_size,
                     plan.counts_a, plan.displs_a,
                     plan.comm, plan.mpi_mode, stream);
    sync_for_timing(timings, stream);
    add_elapsed(timings, &SolveTimings::mpi_backward, mpi_backward_start);

    const double unpack_backward_start = wall_time();
    for (int r = 0; r < plan.nprocs; ++r) {
        const int sub0 = desc_at(plan.gather_nrd_local, 0, r);
        const int sub1 = desc_at(plan.gather_nrd_local, 1, r);
        const int start0 = desc_at(plan.gather_nrd_start, 0, r);
        const int start1 = desc_at(plan.gather_nrd_start, 1, r);
        launch_unpack(plan.drd, plan.nrd_global[0], plan.nrd_global[1],
                      sub0, sub1, start0, start1, plan.bigbuf_a, plan.bigbuf_a_size,
                      plan.displs_a[r], plan.threads_pack, stream);
    }
    sync_for_timing(timings, stream);
    add_elapsed(timings, &SolveTimings::unpack_backward, unpack_backward_start);

    const double update_start = wall_time();
    pascal_update_kernel<<<plan.blocks_tdma, plan.threads_tdma, 0, stream>>>(
        a_dev, b_dev, c_dev, d_dev, plan.drd, nsys, nrow);
    PASCAL_TDMA_CUDA_CHECK(cudaGetLastError());
    sync_for_timing(timings, stream);
    add_elapsed(timings, &SolveTimings::update_compute, update_start);
    add_elapsed(timings, &SolveTimings::total, total_start);
}

}  // namespace

void solve(PascalTdmaPlan& plan,
           double* a_dev,
           double* b_dev,
           double* c_dev,
           double* d_dev,
           const int nsys,
           const int nrow,
           cudaStream_t stream) {
    solve_impl(plan, a_dev, b_dev, c_dev, d_dev, nsys, nrow, nullptr, stream);
}

void solve_profiled(PascalTdmaPlan& plan,
                    double* a_dev,
                    double* b_dev,
                    double* c_dev,
                    double* d_dev,
                    const int nsys,
                    const int nrow,
                    SolveTimings* timings,
                    cudaStream_t stream) {
    if (timings != nullptr) {
        *timings = SolveTimings{};
    }
    solve_impl(plan, a_dev, b_dev, c_dev, d_dev, nsys, nrow, timings, stream);
}

}  // namespace pascal_tdma
