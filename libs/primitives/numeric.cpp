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
void abs(const dlong N, const memory<T> v, memory<T> absv) {

  if (N<=0) return;

  #pragma omp parallel for
  for(dlong n = 0; n < N; ++n){
    absv[n] = std::abs(v[n]);
  }
}

template void abs(const dlong N, const memory<int> v, memory<int> absv);
template void abs(const dlong N, const memory<long long int> v, memory<long long int> absv);
template void abs(const dlong N, const memory<float> v, memory<float> absv);
template void abs(const dlong N, const memory<double> v, memory<double> absv);

template<typename T>
void set(const dlong N, const T val, memory<T> v) {

  if (N<=0) return;

  #pragma omp parallel for
  for(dlong n = 0; n < N; ++n){
    v[n] = val;
  }
}

template void set(const dlong N, const int val, memory<int> v);
template void set(const dlong N, const long long int val, memory<long long int> v);
template void set(const dlong N, const float val, memory<float> v);
template void set(const dlong N, const double val, memory<double> v);

template<typename T>
void range(const dlong N, const T start, const T step, memory<T> v) {
  if (N<=0) return;

  #pragma omp parallel for
  for(dlong n = 0; n < N; ++n){
    v[n] = start + step*n;
  }
}

template void range(const dlong N, const int start, const int step, memory<int> v);
template void range(const dlong N, const long long int start, const long long int step, memory<long long int> v);
template void range(const dlong N, const float start, const float step, memory<float> v);
template void range(const dlong N, const double start, const double step, memory<double> v);

} //namespace prim

} //namespace libp
