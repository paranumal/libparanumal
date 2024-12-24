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
void runLengthEncode(const dlong N,
                     const memory<T> v,
                     dlong& Ngroups,
                     memory<dlong>& offset) {

  if (N<=0) {
    Ngroups=0;
    return;
  }

  if (N==1) {
    Ngroups = 1;
    if (offset.length() < static_cast<size_t>(Ngroups+1)) offset.malloc(Ngroups+1);
    offset[0] = 0; offset[1] = 1;
    return;
  }

  /*Determine how many groups we have*/
  memory<dlong> diff(N);
  adjacentDifferenceFlag(N, v, diff);
  inclusiveScan(N, diff);

  Ngroups = diff[N-1];

  if (offset.length() < static_cast<size_t>(Ngroups+1)) offset.malloc(Ngroups+1);

  offset[0] = 0;

  /*Get the locations where one group transitions to another*/
  #pragma omp parallel for
  for (dlong n=1; n<N; ++n) {
    if(diff[n] - diff[n-1] > 0) {
      offset[diff[n-1]] = n;
    }
  }

  offset[Ngroups] = N;
}

template
void runLengthEncode(const dlong N,
                     const memory<int> v,
                     dlong& Ngroups,
                     memory<dlong>& offset);

template
void runLengthEncode(const dlong N,
                     const memory<long long int> v,
                     dlong& Ngroups,
                     memory<dlong>& offset);

template
void runLengthEncode(const dlong N,
                     const memory<float> v,
                     dlong& Ngroups,
                     memory<dlong>& offset);

template
void runLengthEncode(const dlong N,
                     const memory<double> v,
                     dlong& Ngroups,
                     memory<dlong>& offset);

template<typename T>
void runLengthEncodeConsecutive(const dlong N,
                                const memory<T> v,
                                dlong Ngroups,
                                memory<dlong> offset) {

  if (N<=0) {
    #pragma omp parallel for
    for (dlong n=0; n<Ngroups+1; ++n) { offset[n] = 0; }
    return;
  }

  offset[0] = 0;

  /*Get the locations where one group transitions to another*/
  #pragma omp parallel for
  for (dlong n=0; n<N+1; ++n) {
    const dlong a = (n==0) ? 0 : v[n-1];
    const dlong b = (n==N) ? Ngroups : v[n];
    if(b - a > 0) {
      for (dlong d=a; d<b;++d ) {
        offset[d+1] = n;
      }
    }
  }
}

template
void runLengthEncodeConsecutive(const dlong N,
                                const memory<int> v,
                                dlong Ngroups,
                                memory<dlong> offset);

template
void runLengthEncodeConsecutive(const dlong N,
                                const memory<long long int> v,
                                dlong Ngroups,
                                memory<dlong> offset);

template
void runLengthEncodeConsecutive(const dlong N,
                                const memory<float> v,
                                dlong Ngroups,
                                memory<dlong> offset);

template
void runLengthEncodeConsecutive(const dlong N,
                                const memory<double> v,
                                dlong Ngroups,
                                memory<dlong> offset);

} //namespace prim

} //namespace libp
