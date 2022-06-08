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

namespace libp {

namespace paradogs {

parCSR TentativeProlongator(const dlong Nf,
                            const dlong Nc,
                            platform_t& platform,
                            comm_t comm,
                            memory<hlong>& FineToCoarse,
                            memory<dfloat>& FineNull,
                            memory<dfloat>& CoarseNull) {
  dlong nnz = Nf;
  memory<nonZero_t> entries(nnz);

  hlong localNf=static_cast<hlong>(Nf);
  hlong NfOffsetL=0, NfOffsetU=0;
  comm.Scan(localNf, NfOffsetU);
  NfOffsetL = NfOffsetU-Nf;

  /* Each entry is the CoarseNull vector entry*/
  #pragma omp parallel for
  for (dlong n=0;n<Nf;++n) {
    entries[n].row = n + NfOffsetL;
    entries[n].col = FineToCoarse[n];
    entries[n].val = FineNull[n];
  }

  parCSR T(Nf, Nc,
           nnz, entries,
           platform, comm);
  entries.free();

  /*Create coarse nullvector*/
  CoarseNull.malloc(T.Ncols);

  /*Init coarse null*/
  #pragma omp parallel for
  for (dlong v=0;v<T.Ncols;++v) CoarseNull[v] = 0.0;

  /*Sum columns of T*/
  //add local nonzeros
  for (dlong n=0;n<T.diag.nnz;++n)
    CoarseNull[T.diag.cols[n]] += T.diag.vals[n] * T.diag.vals[n];

  //add nonlocal nonzeros
  for(dlong n=0; n<T.offd.nnz;++n)
    CoarseNull[T.offd.cols[n]] += T.offd.vals[n] * T.offd.vals[n];

  //add the halo values to their origins
  T.halo.Combine(CoarseNull, 1);

  #pragma omp parallel for
  for (dlong n=0;n<Nc;++n)
    CoarseNull[n] = sqrt(CoarseNull[n]);

  //share the results
  T.halo.Exchange(CoarseNull, 1);

  #pragma omp parallel for
  for (dlong n=0;n<T.diag.nnz;++n)
    T.diag.vals[n] /= CoarseNull[T.diag.cols[n]];

  #pragma omp parallel for
  for (dlong n=0;n<T.offd.nnz;++n)
    T.offd.vals[n] /= CoarseNull[T.offd.cols[n]];

  return T;
}

} //namespace paradogs

} //namespace libp
