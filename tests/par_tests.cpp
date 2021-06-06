#include "parallel/algorithm.hpp"
#include "parallel/geomedian.hpp"

#include <iostream>
#include <array>
#include <numeric>
#include <cassert>
#include <algorithm>

using namespace ser;

using Numeric = double;
struct RealHolder {
    Numeric myVal;
};

int main()
{

    const unsigned S = 1e5;

    MPI_Init(nullptr, nullptr);
    int w, r;
    MPI_Comm_size(MPI_COMM_WORLD, &w);
    MPI_Comm_rank(MPI_COMM_WORLD, &r);

    std::vector<RealHolder> x(S);
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd() + r); //Standard mersenne_twister_engine seeded with rd()
    std::lognormal_distribution<> lognormalDistribution(0.0, 1.0);
    for(auto& v : x) {
        v.myVal = (lognormalDistribution(gen));
    }
    MPI_Barrier(MPI_COMM_WORLD);
    auto t1 = MPI_Wtime();
    par::pcout() << par::find_spatial_median(x.begin(), x.end(), 0.01, MPI_COMM_WORLD, [](auto& v){return v.myVal;}, std::nullopt).value_or(-1.0);
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
