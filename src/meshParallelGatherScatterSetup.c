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

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh.h"

// assume nodes locally sorted by rank then global index
// assume gather and scatter are the same sets
ogs_t *meshParallelGatherScatterSetup(mesh_t *mesh,
                                      dlong Nlocal,
                                      dlong *gatherLocalIds,
                                      hlong *gatherBaseIds,
                                      int   *gatherBaseRanks,
                                      int   *gatherHaloFlags,
                                      int verbose) { 

  int rank;
  rank = mesh->rank;

  ogs_t *ogs = (ogs_t*) calloc(1, sizeof(ogs_t));

  // 0. squeeze out negative globalIds
  hlong *baseIds   = (hlong*) calloc(Nlocal,sizeof(hlong));
  dlong *localIds  = (dlong*) calloc(Nlocal,sizeof(dlong));
  int   *baseRanks = (int*)   calloc(Nlocal,sizeof(int));
  int   *haloFlags = (int*)   calloc(Nlocal,sizeof(int));

  dlong Ntotal = 0;
  dlong Nhalo = 0;
  dlong NnonHalo = 0;
  for (dlong n=0;n<Nlocal;n++) {
    if(gatherBaseIds[n]<0) continue;
    localIds[Ntotal]  = gatherLocalIds[n]; 
    baseIds[Ntotal]   = gatherBaseIds[n];   
    baseRanks[Ntotal] = gatherBaseRanks[n]; 
    haloFlags[Ntotal] = gatherHaloFlags[n];
    Ntotal++;

    Nhalo += (gatherHaloFlags[n] == 1);
    NnonHalo += (gatherHaloFlags[n] != 1);
  }

  // ------------------------------------------------------------
  // 1. count number of unique base nodes on this process
  ogs->NhaloGather = 0;
  ogs->NnonHaloGather = 0;
  ogs->NownedHalo = 0;
  for(dlong n=0;n<Ntotal;++n){
    int test = (n==0) ? 1: (baseIds[n] != baseIds[n-1]);
    ogs->NtotalGather += test;
    if (haloFlags[n]==1) ogs->NhaloGather += test;
    if (haloFlags[n]!=1) ogs->NnonHaloGather += test;
    if ((haloFlags[n]==1)&&(baseRanks[n]==rank)) ogs->NownedHalo +=test;
  }

  ogs->haloGatherBaseIds  = (hlong*) calloc(ogs->NhaloGather, sizeof(hlong)); // offset into sorted list of nodes
  ogs->haloGatherOffsets  = (dlong*) calloc(ogs->NhaloGather+1, sizeof(dlong)); // offset into sorted list of nodes
  ogs->haloGatherLocalIds = (dlong*) calloc(Nhalo, sizeof(dlong));

  ogs->ownedHaloGatherIds = (dlong*) calloc(ogs->NhaloGather, sizeof(dlong));  

  ogs->nonHaloGatherBaseIds  = (hlong*) calloc(ogs->NnonHaloGather, sizeof(hlong)); // offset into sorted list of nodes
  ogs->nonHaloGatherOffsets  = (dlong*) calloc(ogs->NnonHaloGather+1, sizeof(dlong)); // offset into sorted list of nodes
  ogs->nonHaloGatherLocalIds = (dlong*) calloc(NnonHalo, sizeof(dlong));

  dlong haloOffset = ogs->NnonHaloGather;

  // only finds bases
  ogs->NhaloGather = 0;
  ogs->NnonHaloGather = 0;
  Nhalo = 0;
  NnonHalo = 0;
  for(dlong n=0;n<Ntotal;++n){
    int test = (n==0) ? 1: (baseIds[n] != baseIds[n-1]);
    
    // increment unique base counter and record index into shuffled list of nodes
    if(haloFlags[n]==1){
      if(test){
        ogs->haloGatherOffsets[ogs->NhaloGather] = Nhalo;
        ogs->haloGatherBaseIds[ogs->NhaloGather] = baseIds[n]+1;
        
        if (baseRanks[n]==rank) {
          ogs->ownedHaloGatherIds[ogs->NhaloGather] = haloOffset++;
        } else {
          ogs->ownedHaloGatherIds[ogs->NhaloGather] = -1;
        }

        ogs->NhaloGather++;
      }
      ogs->haloGatherLocalIds[Nhalo] = localIds[n];
      Nhalo++;
    } else {
      if(test){
        ogs->nonHaloGatherOffsets[ogs->NnonHaloGather] = NnonHalo;
        ogs->nonHaloGatherBaseIds[ogs->NnonHaloGather] = baseIds[n]+1;
        ogs->NnonHaloGather++;
      }
      ogs->nonHaloGatherLocalIds[NnonHalo] = localIds[n];
      NnonHalo++;
    }
  }
  ogs->haloGatherOffsets[ogs->NhaloGather] = Nhalo;
  ogs->nonHaloGatherOffsets[ogs->NnonHaloGather] = NnonHalo;


  // if there are halo nodes to gather
  if(ogs->NhaloGather){

#if 0
    occa::memory o_gatherTmpPinned = mesh->device.mappedAlloc(ogs->NhaloGather*sizeof(dfloat), NULL);
    ogs->haloGatherTmp = (dfloat*) o_gatherTmpPinned.getMappedPointer(); // (char*) calloc(ogs->NhaloGather*sizeof(dfloat), sizeof(char));
#else
    ogs->haloGatherTmp = (dfloat*) occaHostMallocPinned(mesh->device, ogs->NhaloGather*sizeof(dfloat), NULL, ogs->o_haloGatherTmp);
    ogs->o_haloGatherOffsets  = mesh->device.malloc((ogs->NhaloGather+1)*sizeof(dlong), ogs->haloGatherOffsets);
    ogs->o_haloGatherLocalIds = mesh->device.malloc(Nhalo*sizeof(dlong),                ogs->haloGatherLocalIds);
#endif

    ogs->o_ownedHaloGatherIds = mesh->device.malloc(ogs->NhaloGather*sizeof(dlong), ogs->ownedHaloGatherIds);

    // initiate gslib gather-scatter comm pattern on halo nodes only
    ogs->haloGsh = gsParallelGatherScatterSetup(mesh->comm, ogs->NhaloGather, ogs->haloGatherBaseIds,verbose);
  }

  // if there are non-halo nodes to gather
  if(ogs->NnonHaloGather){
    ogs->o_nonHaloGatherOffsets  = mesh->device.malloc((ogs->NnonHaloGather+1)*sizeof(dlong), ogs->nonHaloGatherOffsets);
    ogs->o_nonHaloGatherLocalIds = mesh->device.malloc(NnonHalo*sizeof(dlong),                ogs->nonHaloGatherLocalIds);
  }

  //number of owned gathered nodes
  ogs->Ngather = ogs->NnonHaloGather+ogs->NownedHalo; 

  // build degree vectors
  ogs->invDegree = (dfloat*) calloc(Nlocal, sizeof(dfloat));
  ogs->gatherInvDegree = (dfloat*) calloc(ogs->Ngather, sizeof(dfloat));
  for(dlong n=0;n<Nlocal;++n) ogs->invDegree[n] = 1;

  ogs->o_invDegree = mesh->device.malloc(Nlocal*sizeof(dfloat), ogs->invDegree);
  ogs->o_gatherInvDegree = mesh->device.malloc(ogs->Ngather*sizeof(dfloat), ogs->gatherInvDegree);
  
  meshParallelGather(mesh, ogs, ogs->o_invDegree, ogs->o_gatherInvDegree);

  if(ogs->Ngather)
    ogs->o_gatherInvDegree.copyTo(ogs->gatherInvDegree);

  meshParallelScatter(mesh, ogs, ogs->o_gatherInvDegree, ogs->o_invDegree);

  ogs->o_invDegree.copyTo(ogs->invDegree);  
  
  for(dlong n=0;n<Nlocal;++n) 
    ogs->invDegree[n] = 1./ogs->invDegree[n];
  
  for(dlong n=0;n<ogs->Ngather;++n)
    ogs->gatherInvDegree[n] = 1./ogs->gatherInvDegree[n];

  if(ogs->Ngather)
    ogs->o_gatherInvDegree.copyFrom(ogs->gatherInvDegree);

  if(Nlocal)
    ogs->o_invDegree.copyFrom(ogs->invDegree);

  free(baseIds  );
  free(localIds );
  free(baseRanks);
  free(haloFlags);

  return ogs;
}
