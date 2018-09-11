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

namespace parAlmond {

void formAggregates(parCSR *A, parCSR *C,
                     hlong* FineToCoarse,
                     hlong* globalAggStarts){

  int rank, size;
  MPI_Comm_rank(A->comm, &rank);
  MPI_Comm_size(A->comm, &size);

  const dlong N   = C->Nrows;
  const dlong M   = C->Ncols;
  const dlong diagNNZ = C->diag->nnz;
  const dlong offdNNZ = C->offd->nnz;

  dfloat *rands = (dfloat *) calloc(M, sizeof(dfloat));
  int   *states = (int *)    calloc(M, sizeof(int));

  dfloat *Tr = (dfloat *) calloc(M, sizeof(dfloat));
  int    *Ts = (int *)    calloc(M, sizeof(int));
  hlong  *Ti = (hlong *)  calloc(M, sizeof(hlong));
  hlong  *Tc = (hlong *)  calloc(M, sizeof(hlong));

  hlong *globalRowStarts = A->globalRowStarts;

  for(dlong i=0; i<N; i++)
    rands[i] = (dfloat) drand48();

  // add the number of non-zeros in each column
  int *colCnt = (int *) calloc(M,sizeof(int));
  for(dlong i=0; i<diagNNZ; i++)
    colCnt[C->diag->cols[i]]++;

  for(dlong i=0; i<offdNNZ; i++)
    colCnt[C->offd->cols[i]]++;

  //gs for total column counts
  ogsGatherScatter(colCnt, ogsInt, ogsAdd, A->ogs);

  //add random pertubation
  for(int i=0;i<N;++i)
    rands[i] += colCnt[i];

  //gs to fill halo region
  ogsGatherScatter(rands, ogsDfloat, ogsAdd, A->ogs);

  hlong done = 0;
  while(!done){
    // first neighbours
    // #pragma omp parallel for
    for(dlong i=0; i<N; i++){

      int smax = states[i];
      dfloat rmax = rands[i];
      hlong imax = i + globalRowStarts[rank];

      if(smax != 1){
        //local entries
        for(dlong jj=C->diag->rowStarts[i];jj<C->diag->rowStarts[i+1];jj++){
          const dlong col = C->diag->cols[jj];
          if (col==i) continue;
          if(customLess(smax, rmax, imax, states[col], rands[col], col + globalRowStarts[rank])){
            smax = states[col];
            rmax = rands[col];
            imax = col + globalRowStarts[rank];
          }
        }
        //nonlocal entries
        for(dlong jj=C->offd->rowStarts[i];jj<C->offd->rowStarts[i+1];jj++){
          const dlong col = C->offd->cols[jj];
          if(customLess(smax, rmax, imax, states[col], rands[col], A->colMap[col])) {
            smax = states[col];
            rmax = rands[col];
            imax = A->colMap[col];
          }
        }
      }
      Ts[i] = smax;
      Tr[i] = rmax;
      Ti[i] = imax;
    }

    //share results
    for (dlong n=N;n<M;n++) {
      Tr[n] = 0.;
      Ts[n] = 0;
      Ti[n] = 0;
    }
    ogsGatherScatter(Tr, ogsDfloat, ogsAdd, A->ogs);
    ogsGatherScatter(Ts, ogsInt,    ogsAdd, A->ogs);
    ogsGatherScatter(Ti, ogsHlong,  ogsAdd, A->ogs);

    // second neighbours
    // #pragma omp parallel for
    for(dlong i=0; i<N; i++){
      int    smax = Ts[i];
      dfloat rmax = Tr[i];
      hlong  imax = Ti[i];

      //local entries
      for(dlong jj=C->diag->rowStarts[i];jj<C->diag->rowStarts[i+1];jj++){
        const dlong col = C->diag->cols[jj];
        if (col==i) continue;
        if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
          smax = Ts[col];
          rmax = Tr[col];
          imax = Ti[col];
        }
      }
      //nonlocal entries
      for(dlong jj=C->offd->rowStarts[i];jj<C->offd->rowStarts[i+1];jj++){
        const dlong col = C->offd->cols[jj];
        if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
          smax = Ts[col];
          rmax = Tr[col];
          imax = Ti[col];
        }
      }

      // if I am the strongest among all the 1 and 2 ring neighbours
      // I am an MIS node
      if((states[i] == 0) && (imax == (i + globalRowStarts[rank])))
        states[i] = 1;

      // if there is an MIS node within distance 2, I am removed
      if((states[i] == 0) && (smax == 1))
        states[i] = -1;
    }

    //share results
    for (dlong n=N;n<M;n++) states[n] = 0;
    ogsGatherScatter(states, ogsInt, ogsAdd, A->ogs);

    // if number of undecided nodes = 0, algorithm terminates
    hlong cnt = std::count(states, states+N, 0);
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

  numAggs = 0;
  // enumerate the coarse nodes/aggregates
  for(dlong i=0; i<N; i++) {
    if(states[i] == 1) {
      FineToCoarse[i] = globalAggStarts[rank] + numAggs++;
    } else {
      FineToCoarse[i] = -1;
    }
  }
  for(dlong i=N; i<M; i++) FineToCoarse[i] = 0;

  //share the initial aggregate flags
  ogsGatherScatter(FineToCoarse, ogsHlong, ogsAdd, A->ogs);

  // form the aggregates
  // #pragma omp parallel for
  for(dlong i=0; i<N; i++){
    int   smax = states[i];
    dfloat rmax = rands[i];
    hlong  imax = i + globalRowStarts[rank];
    hlong  cmax = FineToCoarse[i];

    if(smax != 1){
      //local entries
      for(dlong jj=C->diag->rowStarts[i];jj<C->diag->rowStarts[i+1];jj++){
        const dlong col = C->diag->cols[jj];
        if (col==i) continue;
        if(customLess(smax, rmax, imax, states[col], rands[col], col + globalRowStarts[rank])){
          smax = states[col];
          rmax = rands[col];
          imax = col + globalRowStarts[rank];
          cmax = FineToCoarse[col];
        }
      }
      //nonlocal entries
      for(dlong jj=C->offd->rowStarts[i];jj<C->offd->rowStarts[i+1];jj++){
        const dlong col = C->offd->cols[jj];
        if(customLess(smax, rmax, imax, states[col], rands[col], A->colMap[col])){
          smax = states[col];
          rmax = rands[col];
          imax = A->colMap[col];
          cmax = FineToCoarse[col];
        }
      }
    }
    Ts[i] = smax;
    Tr[i] = rmax;
    Ti[i] = imax;
    Tc[i] = cmax;

    if((states[i] == -1) && (smax == 1) && (cmax > -1))
      FineToCoarse[i] = cmax;
  }

  //share results
  for (dlong n=N;n<M;n++) {
    FineToCoarse[n] = 0;
    Tr[n] = 0.;
    Ts[n] = 0;
    Ti[n] = 0;
    Tc[n] = 0;
  }
  ogsGatherScatter(FineToCoarse, ogsHlong,  ogsAdd, A->ogs);
  ogsGatherScatter(Tr,     ogsDfloat, ogsAdd, A->ogs);
  ogsGatherScatter(Ts,     ogsInt,    ogsAdd, A->ogs);
  ogsGatherScatter(Ti,     ogsHlong,  ogsAdd, A->ogs);
  ogsGatherScatter(Tc,     ogsHlong,  ogsAdd, A->ogs);

  // second neighbours
  // #pragma omp parallel for
  for(dlong i=0; i<N; i++){
    int    smax = Ts[i];
    dfloat rmax = Tr[i];
    hlong  imax = Ti[i];
    hlong  cmax = Tc[i];

    //local entries
    for(dlong jj=C->diag->rowStarts[i];jj<C->diag->rowStarts[i+1];jj++){
      const dlong col = C->diag->cols[jj];
      if (col==i) continue;
      if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
        smax = Ts[col];
        rmax = Tr[col];
        imax = Ti[col];
        cmax = Tc[col];
      }
    }
    //nonlocal entries
    for(dlong jj=C->offd->rowStarts[i];jj<C->offd->rowStarts[i+1];jj++){
      const dlong col = C->offd->cols[jj];
      if(customLess(smax, rmax, imax, Ts[col], Tr[col], Ti[col])){
        smax = Ts[col];
        rmax = Tr[col];
        imax = Ti[col];
        cmax = Tc[col];
      }
    }

    if((states[i] == -1) && (smax == 1) && (cmax > -1))
      FineToCoarse[i] = cmax;
  }

  //share results
  for (dlong n=N;n<M;n++) FineToCoarse[n] = 0;
  ogsGatherScatter(FineToCoarse, ogsHlong,  ogsAdd, A->ogs);

  free(rands);
  free(states);
  free(Tr);
  free(Ts);
  free(Ti);
  free(Tc);

  delete C;
}

} //namespace parAlmond