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
                                      int  *gatherHaloFlags,
                                      int verbose) { 

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  ogs_t *ogs = (ogs_t*) calloc(1, sizeof(ogs_t));

  // 0. squeeze out negative globalIds
  iint *baseIds   = (iint*) calloc(Nlocal,sizeof(iint));
  iint *localIds  = (iint*) calloc(Nlocal,sizeof(iint));
  iint *baseRanks = (iint*) calloc(Nlocal,sizeof(iint));
  int  *haloFlags = (iint*) calloc(Nlocal,sizeof(int));

  iint Ntotal = 0;
  iint Nhalo = 0;
  iint NnonHalo = 0;
  for (iint n=0;n<Nlocal;n++) {
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
  for(iint n=0;n<Ntotal;++n){
    int test = (n==0) ? 1: (baseIds[n] != baseIds[n-1]);
    ogs->NtotalGather += test;
    if (haloFlags[n]==1) ogs->NhaloGather += test;
    if (haloFlags[n]!=1) ogs->NnonHaloGather += test;
    if ((haloFlags[n]==1)&&(baseRanks[n]==rank)) ogs->NownedHalo +=test;
  }

  ogs->haloGatherBaseIds  = (iint*) calloc(ogs->NhaloGather, sizeof(iint)); // offset into sorted list of nodes
  ogs->haloGatherOffsets  = (iint*) calloc(ogs->NhaloGather+1, sizeof(iint)); // offset into sorted list of nodes
  ogs->haloGatherLocalIds = (iint*) calloc(Nhalo, sizeof(iint));

  ogs->ownedHaloGatherIds = (iint*) calloc(ogs->NhaloGather, sizeof(iint));  

  ogs->nonHaloGatherBaseIds  = (iint*) calloc(ogs->NnonHaloGather, sizeof(iint)); // offset into sorted list of nodes
  ogs->nonHaloGatherOffsets  = (iint*) calloc(ogs->NnonHaloGather+1, sizeof(iint)); // offset into sorted list of nodes
  ogs->nonHaloGatherLocalIds = (iint*) calloc(NnonHalo, sizeof(iint));

  iint haloOffset = ogs->NnonHaloGather;

  // only finds bases
  ogs->NhaloGather = 0;
  ogs->NnonHaloGather = 0;
  Nhalo = 0;
  NnonHalo = 0;
  for(iint n=0;n<Ntotal;++n){
    iint test = (n==0) ? 1: (baseIds[n] != baseIds[n-1]);
    
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

    occa::memory o_gatherTmpPinned = mesh->device.mappedAlloc(ogs->NhaloGather*sizeof(dfloat), NULL);
    ogs->haloGatherTmp = (dfloat*) o_gatherTmpPinned.getMappedPointer(); // (char*) calloc(ogs->NhaloGather*sizeof(dfloat), sizeof(char));
    
    ogs->o_haloGatherTmp      = mesh->device.malloc(ogs->NhaloGather*sizeof(dfloat),           ogs->haloGatherTmp);
    ogs->o_haloGatherOffsets  = mesh->device.malloc((ogs->NhaloGather+1)*sizeof(iint), ogs->haloGatherOffsets);
    ogs->o_haloGatherLocalIds = mesh->device.malloc(Nhalo*sizeof(iint),                ogs->haloGatherLocalIds);

    ogs->o_ownedHaloGatherIds = mesh->device.malloc(ogs->NhaloGather*sizeof(iint), ogs->ownedHaloGatherIds);

    // initiate gslib gather-scatter comm pattern on halo nodes only
    ogs->haloGsh = gsParallelGatherScatterSetup(ogs->NhaloGather, ogs->haloGatherBaseIds,verbose);
  }

  // if there are non-halo nodes to gather
  if(ogs->NnonHaloGather){
    ogs->o_nonHaloGatherOffsets  = mesh->device.malloc((ogs->NnonHaloGather+1)*sizeof(iint), ogs->nonHaloGatherOffsets);
    ogs->o_nonHaloGatherLocalIds = mesh->device.malloc(NnonHalo*sizeof(iint),                ogs->nonHaloGatherLocalIds);
  }

  //number of owned gathered nodes
  ogs->Ngather = ogs->NnonHaloGather+ogs->NownedHalo; 

  // build degree vectors
  ogs->invDegree = (dfloat*) calloc(Nlocal, sizeof(dfloat));
  ogs->gatherInvDegree = (dfloat*) calloc(ogs->Ngather, sizeof(dfloat));
  for(iint n=0;n<Nlocal;++n) ogs->invDegree[n] = 1;

  ogs->o_invDegree = mesh->device.malloc(Nlocal*sizeof(dfloat), ogs->invDegree);
  ogs->o_gatherInvDegree = mesh->device.malloc(ogs->Ngather*sizeof(dfloat), ogs->gatherInvDegree);
  
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

  free(baseIds  );
  free(localIds );
  free(baseRanks);
  free(haloFlags);

  return ogs;
}
