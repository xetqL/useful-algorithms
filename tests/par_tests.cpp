#include "parallel/algorithm.hpp"
#include "parallel/geomedian.hpp"
#include "parallel/PSRS.hpp"

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

    const unsigned S = 4e7;

    MPI_Init(nullptr, nullptr);
    int w, r;
    MPI_Comm_size(MPI_COMM_WORLD, &w);
    MPI_Comm_rank(MPI_COMM_WORLD, &r);

    MPI_Aint displacements[1]  = {offsetof(RealHolder, myVal)};
    int block_lengths[1]  = {1};
    MPI_Datatype types[1] = {par::get_mpi_type<Numeric>()};
    MPI_Datatype custom_dt;
    MPI_Type_create_struct(1, block_lengths, displacements, types, &custom_dt);
    MPI_Type_commit(&custom_dt);

    std::vector<RealHolder> x(S), y(S);
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd() + r); //Standard mersenne_twister_engine seeded with rd()
    std::lognormal_distribution<> lognormalDistribution(1.0, 1.0);

    for(auto i = 0; i < S; ++i) {
        x[i].myVal = y[i].myVal = (lognormalDistribution(gen));
    }
    //sort_rs(x, custom_dt, MPI_COMM_WORLD, [](auto& x){return x.myVal;});

    MPI_Barrier(MPI_COMM_WORLD);
    auto t1 = MPI_Wtime();
    auto m = par::find_spatial_median(r, w, x.begin(), x.end(), 0.001, MPI_COMM_WORLD, [](auto& x){return x.myVal;}, std::nullopt);
    MPI_Barrier(MPI_COMM_WORLD);
    auto t2 = MPI_Wtime();
    double med;
    typename decltype(x)::iterator it;

    if(m)
        std::tie(med, it) = m.value();



    par::pcout() << std::is_partitioned(x.begin(), x.end(), [med](auto& x){return x.myVal < med;}) <<
                    " " << med << std::endl;

    par::pcout() << "Time for parallel is " << (t2-t1) << " [s]" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);

    t1 = MPI_Wtime();
    std::vector<RealHolder> all(S * w);
    MPI_Gather(y.data(), y.size(), custom_dt, all.data(), y.size(), custom_dt, 0, MPI_COMM_WORLD);
    if(!r) {
        std::nth_element(all.begin(), all.begin() + (S*w) / 2, all.end(), [](const auto& a, const auto& b){return a.myVal < b.myVal; });
        par::pcout() << all.at(S*w/2).myVal << std::endl;
        t2 = MPI_Wtime();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    par::pcout() << "Time for serial is " << (t2-t1) << " [s]" << std::endl;

    MPI_Finalize();
    return 0;
}
