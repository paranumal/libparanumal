/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#ifndef PARALMOND_AMGSETUP_HPP
#define PARALMOND_AMGSETUP_HPP

#include "parAlmond.hpp"
#include "parAlmond/parAlmondMultigrid.hpp"
#include "parAlmond/parAlmondAMGLevel.hpp"


namespace parAlmond {

class strongGraph_t {
public:
  platform_t& platform;
  MPI_Comm comm;
  dlong Nrows=0;
  dlong Ncols=0;
  dlong nnz=0;

  dlong  *rowStarts=nullptr;
  dlong  *cols=nullptr;

  strongGraph_t(dlong N, dlong M, platform_t& _platform, MPI_Comm _comm):
    platform(_platform), comm(_comm), Nrows(N), Ncols(M) {}
  ~strongGraph_t() {
    if (rowStarts) free(rowStarts);
    if (cols) free(cols);
  }
};

amgLevel *coarsenAmgLevel(amgLevel *level, dfloat *null,
                          StrengthType strtype, dfloat theta,
                          AggType aggtype);

strongGraph_t* strongGraph(parCSR *A, StrengthType type, dfloat theta);

void formAggregates(parCSR *A, strongGraph_t *C,
                     hlong* FineToCoarse,
                     hlong* globalAggStarts);

parCSR *tentativeProlongator(parCSR *A, hlong *FineToCoarse,
                            hlong *globalAggStarts, dfloat *null);

parCSR *smoothProlongator(parCSR *A, parCSR *T);

parCSR *transpose(parCSR *A);

parCSR *SpMM(parCSR *A, parCSR *B);

parCSR *galerkinProd(parCSR *A, parCSR *P);

}

#endif
