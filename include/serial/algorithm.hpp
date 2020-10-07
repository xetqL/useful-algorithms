#pragma once

#include <iterator>
#include <algorithm>
#include <vector>

namespace ser {
template<class InputIt>
void partial_sort(InputIt begin, InputIt end, int stride) {
  auto itp = begin, itn = itp;
  while(itn != end) {
      while(itn != end && std::distance(itp, itn) < stride){
          itn++;
      }
      std::sort(itp, itn);
      itp = itn;
  }
}

template<class InputIt>
std::vector<typename std::iterator_traits<InputIt>::value_type> partial_medians(InputIt begin, InputIt end, int stride)  {
  std::vector<typename std::iterator_traits<InputIt>::value_type> medians; medians.reserve(std::distance(begin, end) / stride);
  auto itp = begin, itn = itp;
  while(itn != end) {
      while(itn != end && std::distance(itp, itn) < stride){
          itn++;
      }
      std::sort(itp, itn);    
      medians.push_back(*(itp+(int) std::distance(itp, itn)/2));
      itp = itn;
  }
  return medians;
}
}