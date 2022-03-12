/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#ifndef OGS_UTILS_HPP
#define OGS_UTILS_HPP

#include "ogs.hpp"

namespace libp {

namespace ogs {

struct parallelNode_t{

  dlong localId;    // local node id
  hlong baseId;     // original global index

  dlong newId;         // new global id
  int sign;

  int rank; //original rank
  int destRank; //destination rank

};

template<typename T>
struct ogsType {
  static constexpr Type get();
};

template<> struct ogsType<float> {
  static constexpr Type get() { return Float; }
};
template<> struct ogsType<double> {
  static constexpr Type get() { return Double; }
};
template<> struct ogsType<int> {
  static constexpr Type get() { return Int32; }
};
template<> struct ogsType<long long int> {
  static constexpr Type get() { return Int64; }
};

//permute an array A, according to the ordering returned by P
// i.e. for all n, A[P(n)] <- A[n]
template<typename T, class Order>
void permute(const dlong N, memory<T> A, Order P) {

  for(dlong n=0;n<N;++n) {
    //get what index A[n] should move to
    dlong pn = P(A[n]);
    while (pn!=n) {
      //swap
      std::swap(A[n], A[pn]);
      pn = P(A[n]);
    }
  }
}

} //namespace ogs

} //namespace libp

#endif
