#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh.h"

// assume nodes locally sorted by rank then global index
// assume gather and scatter are the same sets
ogs_t *meshParallelGatherScatterSetup(mesh_t *mesh,
                                      iint Nlocal,
                                      iint *gatherLocalIds,
                                      iint *gatherBaseIds,
                                      iint *gatherBaseRanks,
                                      int  *gatherHaloFlags) { 

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  ogs_t *ogs = (ogs_t*) calloc(1, sizeof(ogs_t));

  // 0. squeeze out negative globalIds
  iint *baseIds   = (iint*) calloc(Nlocal,sizeof(iint));
  iint *localIds  = (iint*) calloc(Nlocal,sizeof(iint));
  iint *baseRanks = (iint*) calloc(Nlocal,sizeof(iint));
  int  *haloFlags = (iint*) calloc(Nlocal,sizeof(int));

  iint Ntotal = 0;
  for (iint n=0;n<Nlocal;n++) {
    if(gatherBaseIds[n]<0) continue;
    localIds[Ntotal]  = gatherLocalIds[n]; 
    baseIds[Ntotal]   = gatherBaseIds[n];   
    baseRanks[Ntotal] = gatherBaseRanks[n]; 
    haloFlags[Ntotal] = gatherHaloFlags[n];
    Ntotal++;
  }

  // ------------------------------------------------------------
  // 1. count number of unique base nodes on this process
  ogs->NtotalGather = 0; // assumes at least one base node
  for(iint n=0;n<Ntotal;++n){
    int test = (n==0) ? 1: (baseIds[n] != baseIds[n-1]);
    ogs->NtotalGather += test;
  }

  ogs->gatherOffsets = (iint*) calloc(ogs->NtotalGather+1, sizeof(iint)); // offset into sorted list of nodes

  // only finds bases
  ogs->NtotalGather = 0; // reset counter
  for(iint n=0;n<Ntotal;++n){
    iint test = (n==0) ? 1: (baseIds[n] != baseIds[n-1]);
    if(test){
      ogs->gatherOffsets[ogs->NtotalGather++] = n;  // increment unique base counter and record index into shuffled list of ndoes
    }
  }
  ogs->gatherOffsets[ogs->NtotalGather] = Ntotal;

  ogs->gatherTmp = (dfloat*) calloc(ogs->NtotalGather, sizeof(dfloat));

  // allocate buffers on DEVICE
  ogs->o_gatherTmp      = mesh->device.malloc(ogs->NtotalGather*sizeof(dfloat), ogs->gatherTmp);
  ogs->o_gatherOffsets  = mesh->device.malloc((ogs->NtotalGather+1)*sizeof(iint), ogs->gatherOffsets);
  ogs->o_gatherLocalIds = mesh->device.malloc(Ntotal*sizeof(iint), localIds);

  ogs->gatherLocalIds  = localIds;
  ogs->gatherBaseIds   = baseIds;
  ogs->gatherHaloFlags = haloFlags;

  // list of nodes to extract from DEVICE gathered array
  ogs->Nhalo = 0;
  ogs->NownedHalo = 0;
  for(iint n=0;n<ogs->NtotalGather;++n){
    for(iint id=ogs->gatherOffsets[n];id<ogs->gatherOffsets[n+1];++id){
      if(haloFlags[id]){ // if any of these are labelled as halo then mark
        ++ogs->Nhalo;
        if (baseRanks[id]==rank) ++ogs->NownedHalo;
        break;
      }
    }
  }

  //number of owned gathered nodes
  ogs->Ngather = ogs->NtotalGather+ogs->NownedHalo-ogs->Nhalo; 

  // set up gather-scatter of halo nodes
  ogs->haloGsh = NULL;

  //node list is sorted by owners, locally-owned nodes first. So the first part of the 
  //  haloIds are the locally-owned halo nodes
  if(ogs->Nhalo){
    ogs->haloLocalIds  = (iint*) calloc(ogs->Nhalo, sizeof(iint));
    ogs->haloGlobalIds = (iint*) calloc(ogs->Nhalo, sizeof(iint));

    ogs->Nhalo = 0;
    for(iint n=0;n<ogs->NtotalGather;++n){
      for(iint id=ogs->gatherOffsets[n];id<ogs->gatherOffsets[n+1];++id){
        if(haloFlags[id]){
          ogs->haloLocalIds[ogs->Nhalo] = n;
          ogs->haloGlobalIds[ogs->Nhalo] = baseIds[id];
          ++ogs->Nhalo;
          break;
        }
      }
    }

    // list of halo node indices
    ogs->o_haloLocalIds = mesh->device.malloc(ogs->Nhalo*sizeof(iint), ogs->haloLocalIds);
    
    ogs->haloTmp = (dfloat*) calloc(ogs->Nhalo, sizeof(dfloat));
    ogs->o_haloTmp  = mesh->device.malloc(ogs->Nhalo*sizeof(dfloat), ogs->haloTmp);

    // initiate gslib gather-scatter comm pattern
    ogs->haloGsh = gsParallelGatherScatterSetup(ogs->Nhalo, ogs->haloGlobalIds);
  }

  // build degree vectors
  ogs->invDegree = (dfloat*) calloc(Nlocal, sizeof(dfloat));
  ogs->gatherInvDegree = (dfloat*) calloc(Nlocal, sizeof(dfloat));
  for(iint n=0;n<Nlocal;++n) ogs->invDegree[n] = 1;

  ogs->o_invDegree = mesh->device.malloc(Nlocal*sizeof(dfloat), ogs->invDegree);
  ogs->o_gatherInvDegree = mesh->device.malloc(Nlocal*sizeof(dfloat), ogs->gatherInvDegree);
  
  meshParallelGather(mesh, ogs, ogs->o_invDegree, ogs->o_gatherInvDegree);

  ogs->o_gatherInvDegree.copyTo(ogs->gatherInvDegree);

  meshParallelScatter(mesh, ogs, ogs->o_gatherInvDegree, ogs->o_invDegree);

  ogs->o_invDegree.copyTo(ogs->invDegree);  
  
  for(iint n=0;n<Nlocal;++n) 
    ogs->invDegree[n] = 1./ogs->invDegree[n];
  
  for(iint n=0;n<ogs->Ngather;++n)
    ogs->gatherInvDegree[n] = 1./ogs->gatherInvDegree[n];

  ogs->o_gatherInvDegree.copyFrom(ogs->gatherInvDegree);
  ogs->o_invDegree.copyFrom(ogs->invDegree);

  return ogs;
}
