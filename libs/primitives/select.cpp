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
void select(const dlong N,
            const memory<T> v,
            const T& val,
                  memory<dlong> ids) {

  if (N<=0) return;

  /*Get the locations where v[] = val*/
  memory<int> predicate(N);
  #pragma omp parallel for
  for (dlong n=0; n<N; ++n) {
    predicate[n] = (v[n] == val) ? 1 : 0;
  }

  inclusiveScan(N, predicate);

  if (predicate[0]==1) {
    ids[0] = 0;
  }

  /*Get the locations where one group transitions to another*/
  #pragma omp parallel for
  for (dlong n=1; n<N; ++n) {
    if(predicate[n] - predicate[n-1] > 0) {
      ids[predicate[n-1]] = n;
    }
  }
}

template
void select(const dlong N,
            const memory<int> v,
            const int& val,
            memory<dlong> ids);

template
void select(const dlong N,
            const memory<long long int> v,
            const long long int& val,
            memory<dlong> ids);

template
void select(const dlong N,
            const memory<float> v,
            const float& val,
            memory<dlong> ids);

template
void select(const dlong N,
            const memory<double> v,
            const double& val,
            memory<dlong> ids);


template<typename T>
void unique(const dlong N,
            const memory<T> v,
                  dlong& Nunique,
                  memory<T>& v_unique) {

  if (N<=0) {
    Nunique=0;
    return;
  }

  /*Get the locations where v[] changes value*/
  memory<int> predicate(N);
  adjacentDifferenceFlag(N, v, predicate);

  Nunique = count(N, predicate, 1);

  memory<dlong> ids(Nunique);
  select(N, predicate, 1, ids);

  if (v_unique.length() < static_cast<size_t>(Nunique)) v_unique.malloc(Nunique);

  transformGather(Nunique, ids, v, v_unique);
}

template
void unique(const dlong N,
            const memory<int> v,
                  dlong& Nunique,
                  memory<int>& v_unique);

template
void unique(const dlong N,
            const memory<long long int> v,
                  dlong& Nunique,
                  memory<long long int>& v_unique);

template
void unique(const dlong N,
            const memory<float> v,
                  dlong& Nunique,
                  memory<float>& v_unique);

template
void unique(const dlong N,
            const memory<double> v,
                  dlong& Nunique,
                  memory<double>& v_unique);


} //namespace prim

} //namespace libp
