//
// Created by xetql on 6/7/21.
//

#pragma once

#include <mpi.h>

template<class T, class GetVal>
void sort_rs(std::vector<T>& data, MPI_Datatype datatype, MPI_Comm comm, GetVal getVal){

    using value_type = typename std::conditional<std::is_arithmetic_v<T>, T, decltype( std::declval<GetVal>()( std::forward<T&>(std::declval<T>())) )>::type;

    MPI_Datatype mpi_value_type = par::get_mpi_type<value_type>();
    int nprocs, rank;

    auto beg = data.begin();
    auto end = data.end();

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    const auto n_samples = nprocs - 1;

    size_t n_local_el = std::distance(beg, end);

    const auto comp = [&getVal](auto a, auto b){ return getVal(a) < getVal(b); };

    // sort n/p*log(n/p)
    std::sort(beg, end, comp);

    // sample p - 1 elements evenly spaced from my list
    std::vector<value_type> my_samples; my_samples.reserve(n_samples);

    for(auto i = 1; i < nprocs; ++i) {
        my_samples.push_back( getVal(beg[(i*n_local_el)/nprocs]) );
    }
    const auto n_all_samples = nprocs * n_samples; //p(p-1) samples
    std::unique_ptr<value_type []> samples_of_samples(new value_type [n_samples]);
    std::unique_ptr<value_type []> all_samples(new value_type [n_all_samples]);

    MPI_Gather(my_samples.data(), n_samples, mpi_value_type, all_samples.get(), n_samples, datatype, 0, comm);

    if(!rank) {
        std::sort(all_samples.get(), all_samples.get()+n_all_samples);
        for(auto i = 0; i < n_samples; ++i) {
            samples_of_samples[i] = all_samples[i * nprocs];
        }
    }

    MPI_Bcast(samples_of_samples.get(), n_samples, mpi_value_type, 0, comm);

    auto prev_it = beg;
    std::vector<MPI_Request > sreqs(nprocs);
    for(auto i = 0; i < n_samples; ++i) {
        auto next_it = std::lower_bound(prev_it, end, samples_of_samples[i], [&getVal](auto& el, auto val){ return getVal(el) < val; });
        auto count = std::distance(prev_it, next_it);
        MPI_Isend(&(*prev_it), count, datatype, i, 9876, comm, &sreqs.at(i));
        prev_it = next_it;
    }
    MPI_Isend(&(*prev_it), std::distance(prev_it, end), datatype, nprocs-1, 9876, comm, &sreqs.at(nprocs-1));

    MPI_Status status;
    std::vector<int> count(nprocs);

    for(auto i = 0; i < nprocs; ++i) {
        MPI_Probe(i, 9876, comm, &status);
        MPI_Get_count(&status, datatype, &count.at(i));
    }

    auto total_count = std::accumulate(count.begin(), count.end(), 0);
    std::unique_ptr<T[]> buff(new T[total_count]);

    for(auto i = 0; i < nprocs; ++i) {
        auto displs = std::accumulate(count.begin(), count.begin() + i, 0);
        MPI_Recv(buff.get() + displs, count.at(i), datatype, i, 9876, comm, MPI_STATUS_IGNORE);
    }

    MPI_Waitall(sreqs.size(), sreqs.data(), MPI_STATUSES_IGNORE);

    data.clear();
    data.reserve(total_count);
    std::copy(buff.get(), buff.get()+total_count, std::back_inserter(data));
    std::sort(data.begin(), data.end(), comp);

}
