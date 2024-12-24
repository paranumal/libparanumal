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

// Check if compiler supports the OMP v5.0 'scan' feature
#if (__GNUC__ >= 10) || (__clang_major__ >= 11) //gcc-10 or higher, or clang-11 or higher
#define HAS_OMP_SCAN
#endif

template<typename T>
void exclusiveScan(const dlong N, memory<T> v) {

  if (N<=0) return;

  T scan_v{0};

#if defined(HAS_OMP_SCAN) && defined(_OPENMP) && !defined(LIBP_DEBUG) && !defined(__clang__)
  /*This looks totally wrong, but is the only way OpenMP will compile
  it, and suprisingly does the right thing. Without OpenMP enabled, however
  this *is* totally wrong*/
  /* NC: Also clang currently has a bug with this type of reduction and fails to compile */
  #pragma omp parallel for reduction(inscan, +:scan_v)
  for(int n = 0; n < N; ++n){
    v[n] = scan_v;
    #pragma omp scan exclusive(scan_v)
    scan_v += v[n];
  }
#else
  /*Right one for the non-openmp case*/
  for(int n = 0; n < N; ++n){
    const T val = v[n];
    v[n] = scan_v;
    scan_v += val;
  }
#endif
}

template void exclusiveScan(const dlong N, memory<int> v);
template void exclusiveScan(const dlong N, memory<long long int> v);
template void exclusiveScan(const dlong N, memory<float> v);
template void exclusiveScan(const dlong N, memory<double> v);


template<typename T>
void exclusiveScan(const dlong N, const memory<T> v, memory<T> w) {

  if (N<=0) return;

  T scan_v{0};

#ifdef HAS_OMP_SCAN
  #pragma omp parallel for reduction(inscan, +:scan_v)
  for(int n = 0; n < N; ++n){
    w[n] = scan_v;
    #pragma omp scan exclusive(scan_v)
    scan_v += v[n];
  }
#else
  for(int n = 0; n < N; ++n){
    w[n] = scan_v;
    scan_v += v[n];
  }
#endif

}

template void exclusiveScan(const dlong N, const memory<int> v, memory<int> w);
template void exclusiveScan(const dlong N, const memory<long long int> v, memory<long long int> w);
template void exclusiveScan(const dlong N, const memory<float> v, memory<float> w);
template void exclusiveScan(const dlong N, const memory<double> v, memory<double> w);

template<typename T>
void inclusiveScan(const dlong N, memory<T> v) {

  if (N<=0) return;

  T scan_v{0};

#ifdef HAS_OMP_SCAN
  #pragma omp parallel for reduction(inscan, +:scan_v)
  for(int n = 0; n < N; ++n){
    scan_v += v[n];
    #pragma omp scan inclusive(scan_v)
    v[n] = scan_v;
  }
#else
  for(int n = 0; n < N; ++n){
    scan_v += v[n];
    v[n] = scan_v;
  }
#endif
}

template void inclusiveScan(const dlong N, memory<int> v);
template void inclusiveScan(const dlong N, memory<long long int> v);
template void inclusiveScan(const dlong N, memory<float> v);
template void inclusiveScan(const dlong N, memory<double> v);


template<typename T>
void inclusiveScan(const dlong N, const memory<T> v, memory<T> w) {

  if (N<=0) return;

  T scan_v{0};

#ifdef HAS_OMP_SCAN
  #pragma omp parallel for reduction(inscan, +:scan_v)
  for(int n = 0; n < N; ++n){
    scan_v += v[n];
    #pragma omp scan inclusive(scan_v)
    w[n] = scan_v;
  }
#else
  for(int n = 0; n < N; ++n){
    scan_v += v[n];
    w[n] = scan_v;
  }
#endif
}

template void inclusiveScan(const dlong N, const memory<int> v, memory<int> w);
template void inclusiveScan(const dlong N, const memory<long long int> v, memory<long long int> w);
template void inclusiveScan(const dlong N, const memory<float> v, memory<float> w);
template void inclusiveScan(const dlong N, const memory<double> v, memory<double> w);


} //namespace prim

} //namespace libp
