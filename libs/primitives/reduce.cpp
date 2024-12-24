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
#include <limits>

namespace libp {

namespace prim {


template<typename T>
T min(const dlong N, const memory<T> v) {

  T min_v = std::numeric_limits<T>::max();

  if (N<=0) return min_v;

  #pragma omp parallel for reduction(min:min_v)
  for(dlong n = 0; n < N; ++n){
    min_v = std::min(min_v, v[n]);
  }

  return min_v;
}

template int min(const dlong N, const memory<int> v);
template long long int min(const dlong N, const memory<long long int> v);
template float min(const dlong N, const memory<float> v);
template double min(const dlong N, const memory<double> v);

template<typename T>
T max(const dlong N, const memory<T> v) {

  T max_v = -std::numeric_limits<T>::max();

  if (N<=0) return max_v;

  #pragma omp parallel for reduction(max:max_v)
  for(dlong n = 0; n < N; ++n){
    max_v = std::max(max_v, v[n]);
  }

  return max_v;
}

template int max(const dlong N, const memory<int> v);
template long long int max(const dlong N, const memory<long long int> v);
template float max(const dlong N, const memory<float> v);
template double max(const dlong N, const memory<double> v);

template<typename T>
T sum(const dlong N, const memory<T> v) {

  T sum_v = T{0};

  if (N<=0) return sum_v;

  #pragma omp parallel for reduction(+:sum_v)
  for(dlong n = 0; n < N; ++n){
    sum_v += v[n];
  }

  return sum_v;
}

template int sum(const dlong N, const memory<int> v);
template long long int sum(const dlong N, const memory<long long int> v);
template float sum(const dlong N, const memory<float> v);
template double sum(const dlong N, const memory<double> v);

} //namespace prim

} //namespace libp
