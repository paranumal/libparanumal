/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "mesh.hpp"

void mesh_t::ParallelGatherScatterSetup() {

  dlong Ntotal = Nverts*(Nelements+totalHaloPairs);

  int *minRank = (int *) malloc(Ntotal*sizeof(int));
  int *maxRank = (int *) malloc(Ntotal*sizeof(int));

  for (dlong i=0;i<Ntotal;i++) {
    minRank[i] = rank;
    maxRank[i] = rank;
  }

  hlong localChange = 0, gatherChange = 1;

  // keep comparing numbers on positive and negative traces until convergence
  while(gatherChange>0){

    // reset change counter
    localChange = 0;

    // send halo data and recv into extension of buffer
    halo->Exchange(minRank, Nverts, ogs_int);
    halo->Exchange(maxRank, Nverts, ogs_int);

    // compare trace vertices
    for(dlong e=0;e<Nelements;++e){
      for(int n=0;n<Nfaces*NfaceVertices;++n){
        dlong id  = e*Nfaces*NfaceVertices + n;
        dlong idM = VmapM[id];
        dlong idP = VmapP[id];

        int minRankM = minRank[idM];
        int minRankP = minRank[idP];

        int maxRankM = maxRank[idM];
        int maxRankP = maxRank[idP];

        if(minRankP<minRankM){
          localChange=1;
          minRank[idM] = minRankP;
        }

        if(maxRankP>maxRankM){
          localChange=1;
          maxRank[idM] = maxRankP;
        }
      }
    }

    // sum up changes
    MPI_Allreduce(&localChange, &gatherChange, 1, MPI_HLONG, MPI_MAX, comm);
  }

  // count elements that contribute to global C0 gather-scatter
  dlong globalCount = 0;
  dlong localCount = 0;
  for(dlong e=0;e<Nelements;++e){
    int isHalo = 0;
    for(int n=0;n<Nverts;++n){
      dlong id = e*Nverts+n;
      if ((minRank[id]!=rank)||(maxRank[id]!=rank)) {
        isHalo = 1;
        break;
      }
    }
    globalCount += isHalo;
    localCount += 1-isHalo;
  }

  globalGatherElementList = (dlong*) calloc(globalCount, sizeof(dlong));
  localGatherElementList  = (dlong*) calloc(localCount, sizeof(dlong));

  globalCount = 0;
  localCount = 0;

  for(dlong e=0;e<Nelements;++e){
    int isHalo = 0;
    for(int n=0;n<Nverts;++n){
      dlong id = e*Nverts+n;
      if ((minRank[id]!=rank)||(maxRank[id]!=rank)) {
        isHalo = 1;
        break;
      }
    }
    if(isHalo){
      globalGatherElementList[globalCount++] = e;
    } else{
      localGatherElementList[localCount++] = e;
    }
  }
  //printf("local = %d, global = %d\n", localCount, globalCount);

  NglobalGatherElements = globalCount;
  NlocalGatherElements = localCount;

  free(minRank);
  free(maxRank);
}
