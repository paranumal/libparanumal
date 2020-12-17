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

#include "parAlmond.hpp"
#include "parAlmond/parAlmondAMGSetup.hpp"

namespace parAlmond {

//create coarsened problem
amgLevel *coarsenAmgLevel(amgLevel *level, dfloat *null,
                          StrengthType strtype, dfloat theta,
                          AggType aggtype){

  int size;
  MPI_Comm_size(level->A->comm, &size);

  strongGraph_t *C = strongGraph(level->A, strtype, theta);

  hlong *FineToCoarse = (hlong *) malloc(level->A->Ncols*sizeof(hlong));
  hlong *globalAggStarts = (hlong *) calloc(size+1,sizeof(hlong));

  formAggregates(level->A, C, FineToCoarse, globalAggStarts);
  delete C;

  // adjustPartition(FineToCoarse, settings);

  parCSR *P;
  parCSR *T = tentativeProlongator(level->A, FineToCoarse, globalAggStarts, null);
  if (aggtype == SMOOTHED) {
    P = smoothProlongator(level->A, T);
    delete T;
  } else {
    P = T;
  }

  // R = P^T
  parCSR *R = transpose(P);

  level->P = P;
  level->R = R;

  parCSR *Acoarse;
  if (aggtype == SMOOTHED) {
    parCSR *AP = SpMM(level->A, P);
    Acoarse = SpMM(R, AP);
    delete AP;
  } else {
    Acoarse = galerkinProd(level->A, P); //specialize for unsmoothed aggregation
  }

  Acoarse->diagSetup();

  amgLevel *coarseLevel = new amgLevel(Acoarse,level->settings);

  //update the number of columns required for this level
  level->Ncols = (level->Ncols > R->Ncols) ? level->Ncols : R->Ncols;
  // coarseLevel->Ncols = (coarseLevel->Ncols > P->Ncols) ? coarseLevel->Ncols : P->Ncols;

  free(FineToCoarse);
  free(globalAggStarts);

  return coarseLevel;
}

} //namespace parAlmond
