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

static bool customLess(int smax, dfloat rmax, hlong imax, int s, dfloat r, hlong i){

  if(s > smax) return true;
  if(smax > s) return false;

  if(r > rmax) return true;
  if(rmax > r) return false;

  if(i > imax) return true;
  if(i < imax) return false;

  return false;
}

/*****************************************************************************/
//
// Parallel Distance 2 Maximal Independant Set (MIS-2) graph partitioning
//
/*****************************************************************************/

void formAggregates(parCSR *A, strongGraph_t *C,
                    hlong* FineToCoarse,
                    hlong* globalAggStarts){

  int rank, size;
  MPI_Comm_rank(A->comm, &rank);
  MPI_Comm_size(A->comm, &size);

  const dlong N   = C->Nrows;
  const dlong M   = C->Ncols;
  const dlong nnz = C->nnz;

  dfloat *rands = (dfloat *) calloc(M, sizeof(dfloat));
  int   *states = (int *)    calloc(M, sizeof(int));
  hlong *colMap = A->colMap; //mapping from local column ids to global ids

  dfloat *Tr = (dfloat *) calloc(M, sizeof(dfloat));
  int    *Ts = (int *)    calloc(M, sizeof(int));
  hlong  *Ti = (hlong *)  calloc(M, sizeof(hlong));
  hlong  *Tc = (hlong *)  calloc(M, sizeof(hlong));

  for(dlong i=0; i<N; i++)
    rands[i] = (dfloat) drand48();

  // add the number of non-zeros in each column
  int *colCnt = (int *) calloc(M,sizeof(int));
  for(dlong i=0; i<nnz; i++)
    colCnt[C->cols[i]]++;

  //gs for total column counts
  A->halo->Combine(colCnt, 1, ogs_int);

  //add random pertubation
  for(int i=0;i<N;++i)
    rands[i] += colCnt[i];

  free(colCnt);

  //gs to fill halo region
  A->halo->Exchange(rands, 1, ogs_dfloat);

  hlong done = 0;
  while(!done){
    // first neighbours
    // #pragma omp parallel for
    for(dlong i=0; i<N; i++){
      int    smax = states[i];
      dfloat rmax = rands[i];
      hlong  imax = colMap[i];

      if(smax != 1){
        for(dlong jj=C->rowStarts[i];jj<C->rowStarts[i+1];jj++){
          const dlong col = C->cols[jj];
          if (col==i) continue;
          if(customLess(smax, rmax, imax, states[col], rands[col], colMap[col])){
            smax = states[col];
            rmax = rands[col];
            imax = colMap[col];
          }
        }
      }
      Ts[i] = smax;
      Tr[i] = rmax;
      Ti[i] = imax;
    }

    //share results
    A->halo->Exchange(Tr, 1, ogs_dfloat);
    A->halo->Exchange(Ts, 1, ogs_int);
    A->halo->Exchange(Ti, 1, ogs_hlong);

    // second neighbours
    // #pragma omp parallel for
    for(dlong i=0; i<N; i++){
      int    smax = Ts[i];
      dfloat rmax = Tr[i];
      hlong  imax = Ti[i];

      for(dlong jj=C->rowStarts[i];jj<C->rowStarts[i+1];jj++){
        const dlong col = C->cols[jj];
        if (col==i) continue;
        if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
          smax = Ts[col];
          rmax = Tr[col];
          imax = Ti[col];
        }
      }

      // if I am the strongest among all the 1 and 2 ring neighbours
      // I am an MIS node
      if((states[i] == 0) && (imax == colMap[i]))
        states[i] = 1;

      // if there is an MIS node within distance 2, I am removed
      if((states[i] == 0) && (smax == 1))
        states[i] = -1;
    }

    //share results
    A->halo->Exchange(states, 1, ogs_int);

    // if number of undecided nodes = 0, algorithm terminates
    hlong cnt = 0;
    for (dlong n=0;n<N;n++) if (states[n]==0) cnt++;

    MPI_Allreduce(&cnt,&done,1,MPI_HLONG, MPI_SUM,A->comm);
    done = (done == 0) ? 1 : 0;
  }

  dlong numAggs = 0;
  dlong *gNumAggs = (dlong *) calloc(size,sizeof(dlong));

  // count the coarse nodes/aggregates
  for(dlong i=0; i<N; i++)
    if(states[i] == 1) numAggs++;

  MPI_Allgather(&numAggs,1,MPI_DLONG,gNumAggs,1,MPI_DLONG,A->comm);

  globalAggStarts[0] = 0;
  for (int r=0;r<size;r++)
    globalAggStarts[r+1] = globalAggStarts[r] + gNumAggs[r];

  free(gNumAggs);

  numAggs = 0;
  // enumerate the coarse nodes/aggregates
  for(dlong i=0; i<N; i++) {
    if(states[i] == 1) {
      FineToCoarse[i] = globalAggStarts[rank] + numAggs++;
    } else {
      FineToCoarse[i] = -1;
    }
  }

  //share the initial aggregate flags
  A->halo->Exchange(FineToCoarse, 1, ogs_hlong);

  // form the aggregates
  // #pragma omp parallel for
  for(dlong i=0; i<N; i++){
    int   smax  = states[i];
    dfloat rmax = rands[i];
    hlong  imax = colMap[i];
    hlong  cmax = FineToCoarse[i];

    if(smax != 1){
      for(dlong jj=C->rowStarts[i];jj<C->rowStarts[i+1];jj++){
        const dlong col = C->cols[jj];
        if (col==i) continue;
        if(customLess(smax, rmax, imax, states[col], rands[col], colMap[col])){
          smax = states[col];
          rmax = rands[col];
          imax = colMap[col];
          cmax = FineToCoarse[col];
        }
      }
    }
    Ts[i] = smax;
    Tr[i] = rmax;
    Ti[i] = imax;
    Tc[i] = cmax;

    if((FineToCoarse[i] == -1) && (smax == 1) && (cmax > -1))
      FineToCoarse[i] = cmax;
  }

  //share results
  A->halo->Exchange(FineToCoarse, 1, ogs_hlong);
  A->halo->Exchange(Tr,     1, ogs_dfloat);
  A->halo->Exchange(Ts,     1, ogs_int);
  A->halo->Exchange(Ti,     1, ogs_hlong);
  A->halo->Exchange(Tc,     1, ogs_hlong);

  // second neighbours
  // #pragma omp parallel for
  for(dlong i=0; i<N; i++){
    int    smax = Ts[i];
    dfloat rmax = Tr[i];
    hlong  imax = Ti[i];
    hlong  cmax = Tc[i];

    for(dlong jj=C->rowStarts[i];jj<C->rowStarts[i+1];jj++){
      const dlong col = C->cols[jj];
      if (col==i) continue;
      if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
        smax = Ts[col];
        rmax = Tr[col];
        imax = Ti[col];
        cmax = Tc[col];
      }
    }

    if((FineToCoarse[i] == -1) && (smax == 1) && (cmax > -1))
      FineToCoarse[i] = cmax;
  }

  //share results
  A->halo->Exchange(FineToCoarse, 1, ogs_hlong);

  free(rands);
  free(states);
  free(Tr);
  free(Ts);
  free(Ti);
  free(Tc);
}

} //namespace parAlmond
