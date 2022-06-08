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

#include "parAdogs.hpp"
#include "parAdogs/parAdogsGraph.hpp"
#include "parAdogs/parAdogsPartition.hpp"
#include <algorithm>
#include <limits>

#ifdef GLIBCXX_PARALLEL
#include <parallel/algorithm>
using __gnu_parallel::partition;
#else
using std::partition;
#endif

namespace libp {

namespace paradogs {

static dfloat Pivot(memory<dfloat>& A,
                    const dlong left,
                    const dlong right,
                    const hlong k,
                    const dfloat min,
                    const dfloat max,
                    comm_t comm) {
  /*Start with guessing a pivot halfway between min and max*/
  const dfloat pivot = (min+max)/2.0;

  /*Bail out if we're looking at a tiny window*/
  constexpr dfloat TOL = (sizeof(dfloat)==8) ? 1.0e-13 : 1.0E-5;
  if (max-min < TOL) return pivot;

  dfloat* Am = partition(A.ptr()+left, A.ptr()+right, [pivot](const dfloat& a){ return a <= pivot; });

  /*Get how many entries are globally <= pivot*/
  hlong localCnt = Am-A.ptr();
  hlong globalCnt = localCnt;
  comm.Allreduce(globalCnt);

  if (globalCnt==k) return pivot;

  if (k<globalCnt) {
    return Pivot(A, left, localCnt, k, min, pivot, comm);
  } else {
    return Pivot(A, localCnt, right, k, pivot, max, comm);
  }
}

/* Given a distributed vector F in comm, find a pivot value,
   such that there are globally k entries of F which are <= pivot. */
dfloat ParallelPivot(const dlong N, memory<dfloat>& F,
                     const hlong k, comm_t comm) {

  /*Make a copy of input vector*/
  memory<dfloat> A(N);
  
  #pragma omp parallel for
  for (dlong n=0;n<N;++n) {
    A[n] = F[n];
  }

  /*Find global minimum/maximum*/
  dfloat globalMin=std::numeric_limits<dfloat>::max();
  dfloat globalMax=std::numeric_limits<dfloat>::min();
  for (dlong n=0;n<N;++n) {
    globalMax = std::max(A[n], globalMax);
    globalMin = std::min(A[n], globalMin);
  }
  comm.Allreduce(globalMin, Comm::Min);
  comm.Allreduce(globalMax, Comm::Max);

  /*Find pivot point via binary search*/
  dfloat pivot = Pivot(A, 0, N, k, globalMin, globalMax, comm);

  return pivot;
}

} //namespace paradogs

} //namespace libp
