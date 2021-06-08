//
// Created by xetql on 6/7/21.
//

#pragma once
#include <mpi.h>

template<class T, class GetVal>
void sort_bs(std::vector<T>& data, MPI_Datatype datatype, MPI_Comm comm, GetVal getVal){

    using value_type = typename std::conditional<std::is_arithmetic_v<T>, T, decltype( std::declval<GetVal>()( std::forward<T&>(std::declval<T>())) )>::type;

    MPI_Datatype mpi_value_type = par::get_mpi_type<value_type>();
    int nprocs, rank;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    int d = static_cast<int>(std::log(nprocs));

    auto beg = data.begin();
    auto end = data.end();
    auto len = std::distance(beg, end);
    std::sort(beg, len, [&getVal](auto& a, auto& b){return getVal(a) < getVal(b);});

    for(auto i = 1;i < d; ++i){
        auto window_id = d-i;
        for(auto j = i-1; j >= 0; --j){

        }
    }


}
