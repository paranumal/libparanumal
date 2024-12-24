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
void transformGather(const dlong N,
                     const memory<dlong> ids,
                     const memory<T> v,
                           memory<T> w) {

  if (N<=0) return;

  #pragma omp parallel for
  for (dlong n=0; n<N; ++n) {
    w[n] = v[ids[n]];
  }
}

template void transformGather(const dlong N, const memory<dlong> id, const memory<int> v, memory<int> w);
template void transformGather(const dlong N, const memory<dlong> id, const memory<long long int> v, memory<long long int> w);
template void transformGather(const dlong N, const memory<dlong> id, const memory<float> v, memory<float> w);
template void transformGather(const dlong N, const memory<dlong> id, const memory<double> v, memory<double> w);


template<typename T>
void transformScatter(const dlong N,
                      const memory<dlong> ids,
                      const memory<T> v,
                            memory<T> w) {

  if (N<=0) return;

  #pragma omp parallel for
  for (dlong n=0; n<N; ++n) {
    w[ids[n]] = v[n];
  }
}

template void transformScatter(const dlong N, const memory<dlong> id, const memory<int> v, memory<int> w);
template void transformScatter(const dlong N, const memory<dlong> id, const memory<long long int> v, memory<long long int> w);
template void transformScatter(const dlong N, const memory<dlong> id, const memory<float> v, memory<float> w);
template void transformScatter(const dlong N, const memory<dlong> id, const memory<double> v, memory<double> w);

} //namespace prim

} //namespace libp
