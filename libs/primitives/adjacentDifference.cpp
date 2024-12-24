/*

The MIT License (MIT)

Copyright (c) 2017-2023 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "primitives.hpp"

namespace libp {

namespace prim {


template<typename T>
void adjacentDifference(const dlong N, const memory<T> v, memory<T> diff) {

  if (N<=0) return;

  diff[0] = v[0];

  #pragma omp parallel for
  for (dlong n=1; n<N; ++n) {
    const T a = v[n-1];
    const T b = v[n];
    diff[n] = b - a;
  }
}

template void adjacentDifference(const dlong N, const memory<int> v, memory<int> diff);
template void adjacentDifference(const dlong N, const memory<long long int> v, memory<long long int> diff);
template void adjacentDifference(const dlong N, const memory<float> v, memory<float> diff);
template void adjacentDifference(const dlong N, const memory<double> v, memory<double> diff);

template<typename T>
void adjacentDifferenceFlag(const dlong N, const memory<T> v, memory<dlong> flag) {

  if (N<=0) return;

  flag[0] = 1;

  #pragma omp parallel for
  for (dlong n=1; n<N; ++n) {
    const T a = v[n-1];
    const T b = v[n];
    flag[n] = (b - a) ? 1 : 0;
  }
}

template void adjacentDifferenceFlag(const dlong N, const memory<int> v, memory<dlong> diff);
template void adjacentDifferenceFlag(const dlong N, const memory<long long int> v, memory<dlong> diff);
template void adjacentDifferenceFlag(const dlong N, const memory<float> v, memory<dlong> diff);
template void adjacentDifferenceFlag(const dlong N, const memory<double> v, memory<dlong> diff);

} //namespace prim

} //namespace libp
