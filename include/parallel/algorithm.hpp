//
// Created by xetql on 10/7/20.
//

#ifndef USEFUL_ALGORITHMS_ALGORITHM_HPP
#define USEFUL_ALGORITHMS_ALGORITHM_HPP
#include <vector>
#include <mpi.h>
#include <algorithm>
#include <serial/algorithm.hpp>

namespace par {

    static std::ostream null(nullptr);

    std::ostream& pcout() {
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
    template<class Iter>
    typename Iter::value_type find_nth(Iter itp, Iter itn, size_t look_for, MPI_Comm comm) {
        return find_nth(itp, itn, look_for, comm, [](auto& v){return v;});
    }
    template<class Iter, class GetValFunc>
    typename Iter::value_type find_nth(Iter itp, Iter itn, size_t look_for, MPI_Comm comm, GetValFunc getValue) {
        using T = typename Iter::value_type;
        int ws, rk;
        MPI_Comm_size(comm, &ws);
        MPI_Comm_rank(comm, &rk);
        size_t cnt;
        std::array<long, 2> split_sizes {};
        std::array<T, 2> pivot_msg {};
        std::vector<T> unfiltered_medians(2*ws);
        std::vector<T> all_medians(ws);
        constexpr T zero = (T) 0, one = (T) 1;
        T pivot;
        do {
            T median_of_medians;
            auto medians = ser::partial_medians(itp, itn, 5);
            std::nth_element(medians.begin(),  medians.begin() + medians.size() / 2, medians.end());
            if(!medians.empty()) {
                pivot_msg.at(0) = one;
                pivot_msg.at(1) = *(medians.begin() + medians.size() / 2);
            } else {
                pivot_msg.at(0) = zero;
            }
            MPI_Gather(pivot_msg.data(), 2, get_mpi_type<T>(), unfiltered_medians.data(), 2, get_mpi_type<T>(), 0, MPI_COMM_WORLD);
            cnt = 0;
            for(size_t i = 0; i < ws; ++i) {
                all_medians.at(cnt) = unfiltered_medians.at(2*i+1);
                cnt += unfiltered_medians.at(2*i);
            }
            if(!rk) {
                std::nth_element(all_medians.begin(), all_medians.begin() + (cnt / 2), all_medians.begin()+cnt);
                pivot = *(all_medians.begin() + (cnt / 2));
            }

            MPI_Bcast(&pivot, 1, get_mpi_type<T>(), 0, MPI_COMM_WORLD);

            auto li = std::partition(itp, itn, [pivot, getValue](const auto& v){ return getValue(v)  < pivot; });
            auto pi = std::partition(li,  itn, [pivot, getValue](const auto& v){ return getValue(v) == pivot; });

            split_sizes = {std::distance(itp, li), std::distance(li, pi)};

            MPI_Allreduce(MPI_IN_PLACE, &split_sizes, 2, get_mpi_type<decltype(split_sizes)::value_type>(), MPI_SUM, MPI_COMM_WORLD);

            if(look_for < split_sizes[0]) {
                itn = li;
            } else if(look_for >= split_sizes[0] + split_sizes[1]){
                itp = pi;
                look_for = look_for - split_sizes[0] - split_sizes[1];
            }

        } while (!(look_for >= split_sizes[0] && look_for < split_sizes[0] + split_sizes[1]));

        return pivot;
    }
}
#endif //USEFUL_ALGORITHMS_ALGORITHM_HPP
