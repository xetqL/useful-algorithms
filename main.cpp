#include "serial_algorithm.hpp"
#include <iostream>
#include <array>

int main()
{
    std::array<int, 0> x = {};
    partial_sort(x.begin(),x.end(), 5);
    auto u = partial_medians(x.begin(), x.end(), 5);
    std::copy(u.begin(), u.end(), std::ostream_iterator<int>(std::cout, " "));
}
