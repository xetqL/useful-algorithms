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
std::optional<double> find_spatial_median(It begin, It end, double tol, MPI_Comm comm, GetValue getVal, std::optional<double> guess) {
    const double epsilon = 1e-9, half = 0.5, acceptable_range_min = half - tol, acceptable_range_max = half + tol;
    unsigned iteration = 0;
    auto n_local = std::distance(begin, end);
    decltype(n_local) n_total;

    MPI_Allreduce(&n_local, &n_total, 1, par::get_mpi_type< decltype(n_local) >(), MPI_SUM, comm);

    unsigned nt_gt = 0, nt_lt = 0;
    double current_cut;
    unsigned current_n_total;

    do {
        MPI_Allreduce(&n_local, &current_n_total, 1, par::get_mpi_type< decltype(n_local) >(), MPI_SUM, comm);
        double mass_center = std::accumulate(begin, end, 0.0, [getVal](const auto& prev, const auto& next){ return prev + getVal(next);}) / static_cast<double>(current_n_total);
        MPI_Allreduce(MPI_IN_PLACE, &mass_center, 1, par::get_mpi_type<decltype(mass_center)>(), MPI_SUM, comm);
        current_cut = iteration == 0 ? guess.value_or(mass_center) : mass_center;

        auto n_lower_than_it = std::partition(begin, end, [getVal, current_cut](const auto& v){ return getVal(v) < current_cut; });

        unsigned current_local_n_lt = std::distance(begin, n_lower_than_it), current_global_n_lt, current_global_n_gt;
        MPI_Allreduce(&current_local_n_lt, &current_global_n_lt,1, par::get_mpi_type<decltype(current_local_n_lt)>(), MPI_SUM, comm);
        current_global_n_gt = current_n_total - current_global_n_lt;

        const double lt_fraction = static_cast<double>(nt_lt + current_global_n_lt) / static_cast<double>(n_total);
        const double gt_fraction = static_cast<double>(nt_gt + current_global_n_gt) / static_cast<double>(n_total);


        if((lt_fraction - epsilon >= acceptable_range_min && lt_fraction + epsilon <= acceptable_range_max) ||
                (gt_fraction - epsilon >= acceptable_range_min && gt_fraction + epsilon <= acceptable_range_max)) { // this is good we stop
            return current_cut;
        }

        if(lt_fraction > gt_fraction) { // keep gt half
            end = n_lower_than_it;
            nt_gt += current_global_n_gt;
        } else { // keep lt half
            begin = n_lower_than_it;
            nt_lt += current_global_n_lt;
        }

        n_local = std::distance(begin, end);
    } while(iteration++ < std::log(n_total));

    return std::nullopt;
}

}
#endif //USEFUL_ALGORITHMS_GEOMEDIAN_HPP
