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

#include "elliptic.h"

// assume nodes locally sorted by rank then global index
// assume gather and scatter are the same sets
void ellipticParallelGatherScatterSetup(elliptic_t* elliptic){  

  mesh_t *mesh = elliptic->mesh;

  // setup occa gather scatter
  int verbose = elliptic->options.compareArgs("VERBOSE","TRUE") ? 1:0;
  mesh->ogs = meshParallelGatherScatterSetup(mesh,mesh->Np*mesh->Nelements,
                                             mesh->gatherLocalIds,
                                             mesh->gatherBaseIds,
                                             mesh->gatherBaseRanks,
                                             mesh->gatherHaloFlags,
                                             verbose);
  elliptic->o_invDegree = mesh->ogs->o_invDegree;

  // count elements that contribute to global C0 gather-scatter
  dlong globalCount = 0;
  dlong localCount = 0;
  for(dlong e=0;e<mesh->Nelements;++e){
    int isHalo = 0;
    for(int n=0;n<mesh->Np;++n){
      if(mesh->globalHaloFlags[e*mesh->Np+n]>0){
        isHalo = 1;
        break;
      }
    }
    globalCount += isHalo;
    localCount += 1-isHalo;
  }

  dlong *globalGatherElementList = (dlong*) calloc(globalCount, sizeof(dlong));
  dlong *localGatherElementList  = (dlong*) calloc(localCount, sizeof(dlong));

  globalCount = 0;
  localCount = 0;

  for(dlong e=0;e<mesh->Nelements;++e){
    int isHalo = 0;
    for(int n=0;n<mesh->Np;++n){
      if(mesh->globalHaloFlags[e*mesh->Np+n]>0){
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

  elliptic->NglobalGatherElements = globalCount;
  elliptic->NlocalGatherElements = localCount;

  if(globalCount)
    elliptic->o_globalGatherElementList =
      mesh->device.malloc(globalCount*sizeof(dlong), globalGatherElementList);

  if(localCount)
    elliptic->o_localGatherElementList =
      mesh->device.malloc(localCount*sizeof(dlong), localGatherElementList);
}
