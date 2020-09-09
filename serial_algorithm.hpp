#pragma once
#include <iterator>
#include <algorithm>

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
