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

//create coarsened problem
amgLevel coarsenAmgLevel(amgLevel& level, memory<dfloat>& null,
                         StrengthType strtype, dfloat theta,
                         AggType aggtype){

  parCSR& A = level.A;

  int size = A.comm.size();

  strongGraph_t C = strongGraph(A, strtype, theta);

  memory<hlong> FineToCoarse(A.Ncols);
  memory<hlong> globalAggStarts(size+1);

  formAggregates(A, C, FineToCoarse, globalAggStarts);

  // adjustPartition(FineToCoarse, settings);

  parCSR P;
  parCSR T = tentativeProlongator(A, FineToCoarse, globalAggStarts, null);
  if (aggtype == SMOOTHED) {
    P = smoothProlongator(A, T);
  } else {
    P = T;
  }

  // R = P^T
  parCSR R = transpose(P);

  level.P = P;
  level.R = R;

  parCSR Acoarse;
  if (aggtype == SMOOTHED) {
    parCSR AP = SpMM(A, P);
    Acoarse = SpMM(R, AP);
  } else {
    Acoarse = galerkinProd(A, P); //specialize for unsmoothed aggregation
  }

  Acoarse.diagSetup();

  amgLevel coarseLevel(Acoarse,level.settings);

  //update the number of columns required for this level
  level.Ncols = std::max(level.Ncols, std::max(A.Ncols, R.Ncols));
  // coarseLevel->Ncols = (coarseLevel->Ncols > P->Ncols) ? coarseLevel->Ncols : P->Ncols;

  return coarseLevel;
}

} //namespace parAlmond

} //namespace libp
