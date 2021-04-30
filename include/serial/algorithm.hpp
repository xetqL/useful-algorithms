#pragma once

#include <iterator>
#include <algorithm>
#include <vector>

namespace ser {

template<class ForwardIt, class OutputIt, class BinaryOp>
void combine(ForwardIt beginA, ForwardIt endA, ForwardIt beginB, OutputIt out, BinaryOp op) {
    for(; beginA != endA; beginA++, beginB++) {
        *out = op(*beginA, *beginB);
    }
}

template<class InputIt, class Comp>
void stride_sort(InputIt begin, InputIt end, long stride, Comp comp) {
  auto itp = begin, itn = itp;
  while(itn != end) {
      std::advance(itn, std::min(stride, std::distance(itn, end)));
      std::sort(itp, itn, comp);
      itp = itn;
  }
}

template<class InputIt>
void stride_sort(InputIt begin, InputIt end, long stride){
    stride_sort(begin, end, stride, std::less<>());
}

template<class InputIt, class Comp>
std::vector<typename std::iterator_traits<InputIt>::value_type> stride_medians(InputIt begin, InputIt end, long stride, Comp comp)  {
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
std::vector<typename std::iterator_traits<InputIt>::value_type> stride_medians(InputIt begin, InputIt end, long stride){
    return stride_medians(begin, end, stride, std::less<>());
}
}