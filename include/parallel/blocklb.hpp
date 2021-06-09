//
// Created by xetql on 6/9/21.
//
#include <mpi.h>
#include <algorithm>
#include <vector>
#include <cassert>

#pragma once

template<class T>
void blocklb(int rank, int nprocs, std::vector<T>& data, MPI_Datatype datatype, MPI_Comm comm) {
    auto n_local = data.size();
    decltype(n_local) n_up_to_mine = 0;
    std::vector<size_t> data_ids(n_local);
    decltype(n_local) n_total;

    MPI_Allreduce(&n_local, &n_total,   1, par::get_mpi_type<decltype(n_local)>(), MPI_SUM, comm);
    MPI_Exscan(&n_local, &n_up_to_mine, 1, par::get_mpi_type<decltype(n_local)>(), MPI_SUM, comm);

    int p_over_n = static_cast<double>(nprocs) / static_cast<double>(n_total),
        n_over_p = static_cast<double>(n_total) / static_cast<double>(nprocs);

    int start_id_PE = n_up_to_mine * p_over_n;
    std::vector<int> counts(nprocs), displs(nprocs), rcounts(nprocs), rdispls(nprocs);
    long current_index_lo = start_id_PE * n_over_p - n_up_to_mine;
    long current_index_hi = (start_id_PE+1) * n_over_p - n_up_to_mine;
    for(auto PE = start_id_PE; PE < nprocs &&
        (current_index_lo < static_cast<long>(n_local) || current_index_hi < static_cast<long>(n_local)); ++PE) {

        if(current_index_lo < 0 && current_index_hi <= n_local){
            counts.at(PE) = current_index_hi;
        } else if(current_index_hi > n_local && current_index_lo < n_local) {
            counts.at(PE) = n_local - std::max(static_cast<decltype(current_index_lo)>(0), current_index_lo);
        } else if( current_index_lo >= 0 && current_index_hi >= 0 )
            counts.at(PE) = n_over_p;

        current_index_lo += n_over_p;
        current_index_hi += n_over_p;

    }
    std::exclusive_scan(counts.begin(), counts.end(), displs.begin(), 0, std::plus{});
    MPI_Alltoall(counts.data(), 1, MPI_INT, rcounts.data(), 1, MPI_INT, comm);
    std::exclusive_scan(rcounts.begin(), rcounts.end(), rdispls.begin(), 0, std::plus{});
    std::vector<T> recv_buf(std::accumulate(rcounts.begin(), rcounts.end(), 0));
    MPI_Alltoallv(data.data(), counts.data(), displs.data(), datatype,
              recv_buf.data(), rcounts.data(), rdispls.data(), datatype, comm);

    data = std::move(recv_buf);
}