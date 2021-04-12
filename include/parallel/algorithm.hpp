//
// Created by xetql on 10/7/20.
//

#ifndef USEFUL_ALGORITHMS_ALGORITHM_HPP
#define USEFUL_ALGORITHMS_ALGORITHM_HPP
#include <vector>
#include <mpi.h>
#include <algorithm>
#include <serial/algorithm.hpp>
#include <random>
#include <serial/type.hpp>
#include <serial/format.hpp>

namespace par {

    static std::ostream null(nullptr);

    inline std::ostream& pcout() {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(rank == 0) return std::cout;
        else return null;
    }

    template<class T>
    constexpr MPI_Datatype get_mpi_type() {
        if constexpr (std::is_same<T, float>::value) return MPI_FLOAT;
        if constexpr (std::is_same<T, double>::value) return MPI_DOUBLE;
        if constexpr (std::is_same<T, int>::value) return MPI_INT;
        if constexpr (std::is_same<T, unsigned int>::value) return MPI_UNSIGNED;
        if constexpr (std::is_same<T, long>::value) return MPI_LONG;
        if constexpr (std::is_same<T, long int>::value) return MPI_LONG_INT;
        if constexpr (std::is_same<T, long double>::value) return MPI_LONG_DOUBLE;
        if constexpr (std::is_same<T, long long>::value) return MPI_LONG_LONG;
        if constexpr (std::is_same<T, long long int>::value) return MPI_LONG_LONG_INT;
        if constexpr (std::is_same<T, unsigned long>::value) return MPI_UNSIGNED_LONG;
        if constexpr (std::is_same<T, unsigned long long>::value) return MPI_UNSIGNED_LONG_LONG;
        if constexpr (std::is_same<T, short>::value) return MPI_SHORT;
        if constexpr (std::is_same<T, short int>::value) return MPI_SHORT_INT;
        if constexpr (std::is_same<T, char>::value) return MPI_CHAR;
        return MPI_DATATYPE_NULL;
    }

namespace {
    template<class Iter, class LesserThanComp, class EqComp>
    auto _find_nth(Iter itp, Iter itn, size_t look_for, MPI_Datatype datatype, MPI_Comm comm, LesserThanComp lt, EqComp eq) -> typename Iter::value_type {
        using T = typename Iter::value_type;
        int ws, rk;
        MPI_Comm_size(comm, &ws);
        MPI_Comm_rank(comm, &rk);
        unsigned long split_size[3];
        std::vector<long> pivot_msgs(ws);
        srand(rk + time(nullptr));
        T* pivot = (T*) new char[sizeof(T)];
        T* all_medians = (T*) new char[sizeof(T)*ws];
        do {
            long N = std::distance(itp, itn);
            MPI_Gather(&N, 1, get_mpi_type<long>(), pivot_msgs.data(), 1, get_mpi_type<long>(), 0, comm);

            if (N) {
                T local_pivot = *(itp + (rand() % N));
                MPI_Send(&local_pivot, 1, datatype, 0, 777, comm);
            }

            if(!rk) {
                auto pivot_count = std::count_if(pivot_msgs.begin(), pivot_msgs.end(), [](auto v){return v;});
                for(auto i = 0; i < pivot_count; i++) {
                    MPI_Recv((all_medians+i), 1, datatype, MPI_ANY_SOURCE, 777, comm, MPI_STATUSES_IGNORE);
                }
                std::nth_element(all_medians, all_medians + (pivot_count / 2), all_medians + pivot_count, lt);
                *pivot = all_medians[pivot_count / 2];
            }

            MPI_Bcast(pivot, 1, datatype, 0, comm);

            auto ilt = std::partition(itp, itn, [pivot, lt](const auto& v){ return lt(v, *pivot); });
            auto ieq = std::partition(ilt, itn,  [pivot, eq](const auto& v) { return eq(v, *pivot); });

            split_size[0] = std::distance(itp, ilt);
            split_size[1] = std::distance(ilt, ieq);

            MPI_Allreduce(MPI_IN_PLACE, &split_size, 2, get_mpi_type<unsigned long>(), MPI_SUM, comm);

            auto nlt = split_size[0];
            auto neq = split_size[1];

            if(look_for < nlt) {
                itn = ilt;
            } else if(look_for < (nlt + neq) ) {
                auto retval = *pivot;

                delete[] all_medians;
                delete[] pivot;

                return retval;
            } else {
                itp = ieq;
                look_for = look_for-(nlt+neq);
            }

        } while (true);


    }
}

    template<class Iter, class LtComp, class EqComp>
    typename Iter::value_type find_nth(Iter itp, Iter itn, size_t look_for, MPI_Datatype datatype, MPI_Comm comm, LtComp lt, EqComp eq) {
        int ws;
        MPI_Comm_size(comm, &ws);
        if (ws > 1){
            return _find_nth(itp, itn, look_for, datatype, comm, lt, eq);
        } else {
            std::nth_element(itp, itp + look_for, itn, lt);
            return *(itp+look_for);
        }
    }

    template<class Iter>
    typename Iter::value_type find_nth(Iter itp, Iter itn, size_t look_for, MPI_Comm comm) {
        return find_nth(itp, itn, look_for, get_mpi_type<typename Iter::value_type>(), comm, std::less<>{}, std::equal_to<>{});
    }
}
#endif //USEFUL_ALGORITHMS_ALGORITHM_HPP
