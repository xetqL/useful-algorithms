#include "parallel/algorithm.hpp"

#include <iostream>
#include <array>
#include <numeric>
#include <cassert>
#include <algorithm>

using namespace ser;
int main()
{
    using Numeric = long;

    const Numeric S = 10'000'005;

    MPI_Init(nullptr, nullptr);
    int w, r;
    MPI_Comm_size(MPI_COMM_WORLD, &w);
    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    std::vector<Numeric> x(S);
    std::iota(x.begin(), x.end(), r * S);

    MPI_Barrier(MPI_COMM_WORLD);
    auto t1 = MPI_Wtime();
    par::pcout() << par::find_nth(std::begin(x), std::end(x), (x.size() * w)/2, MPI_COMM_WORLD) << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    auto t2 = MPI_Wtime();
    par::pcout() << "Time for parallel is " << (t2-t1) << " [s]" << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();
    std::vector<Numeric> all(S * w);
    MPI_Gather(x.data(), x.size(), par::get_mpi_type<Numeric>(), all.data(), x.size(), par::get_mpi_type<Numeric>(), 0, MPI_COMM_WORLD);
    if(!r){
        std::nth_element(all.begin(), all.begin()+(w * S / 2), all.end());
        std::cout << *(all.begin()+(w * S / 2)) << std::endl;
        t2 = MPI_Wtime();
    }
    par::pcout() << "Time for serial is " << (t2-t1) << " [s]" << std::endl;

    MPI_Finalize();
}
