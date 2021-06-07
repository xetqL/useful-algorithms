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
    const unsigned S = 100'000;

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
        x[i].myVal  = (lognormalDistribution(gen));
    }

    MPI_Barrier(MPI_COMM_WORLD);
    auto t1 = MPI_Wtime();
    sort_rs(x, custom_dt, MPI_COMM_WORLD, [](auto& v){return v.myVal;});
    MPI_Barrier(MPI_COMM_WORLD);
    auto t2 = MPI_Wtime();

    par::pcout() << "Time for parallel is " << (t2-t1) << " [s]" << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    for(auto i = 0; i < w;++i) {
        for(auto & j : x) {
            if(r == i)
                std::cout << j.myVal << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);

    }

    t1 = MPI_Wtime();
    std::vector<RealHolder> all(S * w);
    //MPI_Gather(x.data(), x.size(), custom_dt, all.data(), x.size(), custom_dt, 0, MPI_COMM_WORLD);
    MPI_Request req;
    MPI_Isend(y.data(), y.size(), custom_dt, 0, 8888, MPI_COMM_WORLD, &req);
    if(!r) {
        int recvd = 0;
        for(auto i = 0; i < w; ++i){
            MPI_Status status;
            MPI_Probe(i , 8888, MPI_COMM_WORLD, &status);
            int count;
            MPI_Get_count(&status, custom_dt, &count);
            MPI_Recv(all.data() + recvd, count, custom_dt, i, 8888, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            recvd += count;
        }
        std::nth_element(all.begin(), all.begin() + ((S*w) / 2), all.end(), [](const auto& a, const auto& b){return a.myVal < b.myVal; });
        t2 = MPI_Wtime();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    par::pcout() << "Time for serial is " << (t2-t1) << " [s]" << std::endl;

    MPI_Finalize();
    return 0;
}
