#include "parallel/algorithm.hpp"

#include <iostream>
#include <array>
#include <numeric>
#include <cassert>
#include <algorithm>

using namespace ser;

struct IntHolder {
    int myVal;
};

int main()
{
    using Numeric = int;

    const Numeric S = 1e4;

    MPI_Init(nullptr, nullptr);
    int w, r;
    MPI_Comm_size(MPI_COMM_WORLD, &w);
    MPI_Comm_rank(MPI_COMM_WORLD, &r);

    MPI_Aint displacements[1]  = {offsetof(IntHolder, myVal)};
    int block_lengths[1]  = {1};
    MPI_Datatype types[1] = {MPI_INT};
    MPI_Datatype custom_dt;
    MPI_Type_create_struct(1, block_lengths, displacements, types, &custom_dt);
    MPI_Type_commit(&custom_dt);

    std::vector<IntHolder> x(S);
    srand(0);
    for(auto& v : x) {
        v.myVal += (int) rand() % 10000;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    auto t1 = MPI_Wtime();
    par::pcout() << par::find_nth(std::begin(x), std::end(x), (x.size() * w)/2, custom_dt, MPI_COMM_WORLD,
                                  [](auto& lhs, auto& rhs){return lhs.myVal < rhs.myVal; },
                                  [](auto& lhs, auto& rhs){return lhs.myVal == rhs.myVal; }
                                  ).myVal << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    auto t2 = MPI_Wtime();
    par::pcout() << "Time for parallel is " << (t2-t1) << " [s]" << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();
    std::vector<Numeric> all(S * w);
    MPI_Gather(x.data(), x.size(), par::get_mpi_type<Numeric>(), all.data(), x.size(), par::get_mpi_type<Numeric>(), 0, MPI_COMM_WORLD);
    if(!r) {
        std::nth_element(all.begin(), all.begin()+(w * S / 2), all.end());
        std::cout << *(all.begin()+(w * S / 2)) << std::endl;
        t2 = MPI_Wtime();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    par::pcout() << "Time for serial is " << (t2-t1) << " [s]" << std::endl;

    MPI_Finalize();
    return 0;
}
