/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#include "parAlmond.hpp"
#include "parAlmond/parAlmondAMGSetup.hpp"

namespace libp {

namespace parAlmond {

parCSR tentativeProlongator(parCSR& A, memory<hlong> FineToCoarse,
                            memory<hlong> globalAggStarts, memory<dfloat> null){

  int rank = A.comm.rank();
  // int size = A.comm.size();

  const dlong NCoarse = static_cast<dlong>(globalAggStarts[rank+1]-globalAggStarts[rank]); //local num agg

  parCOO cooP(A.platform, A.comm);

  //copy global partition
  cooP.globalRowStarts = A.globalRowStarts;
  cooP.globalColStarts = globalAggStarts;

  const hlong globalRowOffset = A.globalRowStarts[rank];

  cooP.nnz = A.Nrows;
  cooP.entries.malloc(cooP.nnz);

  for(dlong n=0; n<cooP.nnz; n++) {
    cooP.entries[n].row = n + globalRowOffset;
    cooP.entries[n].col = FineToCoarse[n];
    cooP.entries[n].val = null[n];
  }

  //build P from coo matrix
  parCSR P(cooP);

  // normalize the columns of P and fill null with coarse null vector

  //check size. If this ever triggers, we'll have to implement a re-alloc of null
  LIBP_ABORT("Size of Coarse nullvector is too large, need to re-alloc",
             P.Ncols > A.Ncols);

  //set coarse null to 0
  for(dlong i=0; i<P.Ncols; i++) null[i] = 0.0;

  //add local nonzeros
  for(dlong i=0; i<P.diag.nnz; i++)
    null[P.diag.cols[i]] += P.diag.vals[i] * P.diag.vals[i];

  //add nonlocal nonzeros
  for(dlong i=0; i<P.offd.nnz; i++)
    null[P.offd.cols[i]] += P.offd.vals[i] * P.offd.vals[i];

  //add the halo values to their origins
  P.halo.Combine(null, 1);

  for(dlong i=0; i<NCoarse; i++)
    null[i] = sqrt(null[i]);

  //share the results
  P.halo.Exchange(null, 1);

  for(dlong i=0; i<P.diag.nnz; i++)
    P.diag.vals[i] /= null[P.diag.cols[i]];
  for(dlong i=0; i<P.offd.nnz; i++)
    P.offd.vals[i] /= null[P.offd.cols[i]];

  return P;
}

} //namespace parAlmond

} //namespace libp
