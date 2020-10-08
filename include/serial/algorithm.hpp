#pragma once

#include <iterator>
#include <algorithm>
#include <vector>

namespace ser {

template<class T>
struct identity_of{
    T operator()(T&& v){
        return std::forward<T>(v);
    }
    T operator()(T& v){
        return std::forward<T>(v);
    }
};
template<class InputIt, class Comp>
void partial_sort(InputIt begin, InputIt end, long stride, Comp comp) {
  auto itp = begin, itn = itp;
  while(itn != end) {
      std::advance(itn, std::min(stride, std::distance(itn, end)));
      std::sort(itp, itn, comp);
      itp = itn;
  }
}
template<class InputIt>
void partial_sort(InputIt begin, InputIt end, long stride){
    partial_sort(begin, end, stride, std::less<>());
}
template<class InputIt, class Comp>
std::vector<typename std::iterator_traits<InputIt>::value_type> partial_medians(InputIt begin, InputIt end, long stride, Comp comp)  {
  std::vector<typename std::iterator_traits<InputIt>::value_type> medians; medians.reserve(std::distance(begin, end) / stride);
  auto itp = begin, itn = itp;
  while(itn != end) {
      std::advance(itn, std::min(stride, std::distance(itn, end)));
      std::sort(itp, itn, comp);
      medians.push_back(*(itp+(std::distance(itp, itn)/2l)));
      itp = itn;
  }
  return medians;
}
template<class InputIt>
std::vector<typename std::iterator_traits<InputIt>::value_type> partial_medians(InputIt begin, InputIt end, long stride){
    return partial_medians(begin, end, stride, std::less<>());
}
}