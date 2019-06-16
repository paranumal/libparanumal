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

  dlong Ntotal = Np*Nelements;

  int verbose = 0;
  ogs = ogsSetup(Ntotal, globalIds, comm, verbose, device);

  //use the gs to find what nodes are local to this rank
  int *minRank = (int *) calloc(Ntotal,sizeof(int));
  int *maxRank = (int *) calloc(Ntotal,sizeof(int));
  for (dlong i=0;i<Ntotal;i++) {
    minRank[i] = rank;
    maxRank[i] = rank;
  }

  ogsGatherScatter(minRank, ogsInt, ogsMin, ogs); //minRank[n] contains the smallest rank taking part in the gather of node n
  ogsGatherScatter(maxRank, ogsInt, ogsMax, ogs); //maxRank[n] contains the largest rank taking part in the gather of node n

  // count elements that contribute to global C0 gather-scatter
  dlong globalCount = 0;
  dlong localCount = 0;
  for(dlong e=0;e<Nelements;++e){
    int isHalo = 0;
    for(int n=0;n<Np;++n){
      dlong id = e*Np+n;
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
    for(int n=0;n<Np;++n){
      dlong id = e*Np+n;
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
}
