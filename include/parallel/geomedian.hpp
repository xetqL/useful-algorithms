//
// Created by xetql on 6/4/21.
//

#ifndef USEFUL_ALGORITHMS_GEOMEDIAN_HPP
#define USEFUL_ALGORITHMS_GEOMEDIAN_HPP
#include <mpi.h>
#include <iterator>
#include "algorithm.hpp"
namespace par {
// can be downcast l8r

template<class It, class GetValue>
std::optional<double> find_spatial_median(It begin, It end, double tol, MPI_Comm comm, GetValue getVal, std::optional<double> guess){
    int rank, nprocs;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);
    auto n_local = std::distance(begin, end);
    decltype(n_local) n_total, n_avg, n_max;

    std::sort(begin, end, [&getVal](auto& a, auto& b){ return getVal(a) < getVal(b); });

    MPI_Allreduce(&n_local, &n_total, 1, par::get_mpi_type< decltype(n_local) >(), MPI_SUM, comm);
    MPI_Allreduce(&n_local, &n_max,   1, par::get_mpi_type< decltype(n_local) >(), MPI_MAX, comm);
    n_avg = n_total / nprocs;
    double imb = static_cast<double>(n_max) / static_cast<double>(n_avg);

    double mass_center = std::accumulate(begin, end, 0.0, [&getVal](auto& prev, auto& next){ return prev + getVal(next);});
    mass_center /= static_cast<double>(n_total);
    MPI_Allreduce(MPI_IN_PLACE, &mass_center, 1, par::get_mpi_type<decltype(mass_center)>(), MPI_SUM, comm);

    // check if actual median lies outside of tolerance zone
    unsigned n_lower_than, n_greater_than;
    if(imb > tol) {
        double current_cut = guess.value_or(mass_center);
        while(true) {
            std::cout << current_cut << std::endl;

            MPI_Bcast(&current_cut, 1, MPI_DOUBLE, 0, comm);
            auto n_lower_than_it = std::lower_bound(begin, end, current_cut, [&getVal](auto& v, auto& value){ return getVal(v) < value; });
            unsigned currennt_n_lower_than;
            MPI_Allreduce(MPI_IN_PLACE, &currennt_n_lower_than, 1, par::get_mpi_type<decltype(n_lower_than)>(), MPI_SUM, comm);
            n_greater_than = n_total - n_lower_than;
            imb = std::fabs(n_lower_than - n_greater_than) / n_total;
            if(imb <= tol) return current_cut;

            if (n_lower_than > n_greater_than) {
                end = n_lower_than_it;
            } else if (n_greater_than > n_lower_than) {
                begin = n_lower_than_it;
            }

            mass_center = std::accumulate(begin, end, 0.0, [&getVal](auto& prev, auto& next){ return prev + getVal(next);});
            mass_center /= static_cast<double>(n_total);
            MPI_Allreduce(MPI_IN_PLACE, &mass_center, 1, par::get_mpi_type<decltype(mass_center)>(), MPI_SUM, comm);
            if(std::fabs(current_cut - mass_center) <= 1e-12) return std::nullopt;
            current_cut = mass_center;
        }
    } else {
        return std::nullopt;
    }
}

}
#endif //USEFUL_ALGORITHMS_GEOMEDIAN_HPP
