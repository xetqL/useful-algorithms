#include "serial/algorithm.hpp"

#include <iostream>
#include <array>
#include <numeric>
#include <cassert>
#include <algorithm>

using namespace ser;
int main()
{
    std::vector<int> x = {3,1,4,2,5, 4,4,4,3,1 , 5,1,0,0,5};
    stride_sort(x.begin(),x.end(), 5);
    auto u = stride_medians(x.begin(), x.end(), 5);
    auto res = {3, 4, 1};
    assert( std::equal(u.begin(), u.end(), res.begin()) );
}
