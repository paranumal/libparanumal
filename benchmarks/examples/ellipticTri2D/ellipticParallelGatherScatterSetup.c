#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include <mpi.h>

#include "mesh.h"

// assume nodes locally sorted by rank then global index
// assume gather and scatter are the same sets
void ellipticParallelGatherScatterSetup(mesh_t *mesh,    // provides DEVICE
                                        iint Nlocal,     // number of local nodes
                                        iint Nbytes,     // number of bytes per node
                                        iint *gatherLocalIds,  // local index of nodes
                                        iint *gatherBaseIds,   // global index of their base nodes
                                        iint *gatherHaloFlags,
                                        ogs_t **halo,
                                        ogs_t **nonHalo){   // 1 for halo node, 0 for not
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#if 1
  // already done in first PGS
  // ------------------------------------------------------------  
  // 0. propagate halo flags uniformly using a disposable gs instance

  void *allGsh = gsParallelGatherScatterSetup(Nlocal, gatherBaseIds);

  // compute max halo flag using global numbering
  gsParallelGatherScatter(allGsh, gatherHaloFlags, "int", "max"); // should use iint

  // tidy up
  gsParallelGatherScatterDestroy(allGsh);
  if(rank==0)
    printf("finished temporary GS setup\n");
#endif
  // initialize gather structs
  *halo    = (ogs_t*) calloc(1, sizeof(ogs_t));
  *nonHalo = (ogs_t*) calloc(1, sizeof(ogs_t));
  
  // ------------------------------------------------------------
  // 1. count number of unique base nodes on this process 
  (*halo)->Ngather = 0;
  (*nonHalo)->Ngather = 0;

  iint nHalo = 0;
  iint nNonHalo = 0;
  
  for(iint n=0;n<Nlocal;++n){
    iint test = (n==0) ? 1: (gatherBaseIds[n] != gatherBaseIds[n-1]);
    if(gatherHaloFlags[n]==1){
      (*halo)->Ngather += test;
      ++nHalo;
    }
  }


  for(iint n=0;n<Nlocal;++n){
    iint test = (n==0) ? 1: (gatherBaseIds[n] != gatherBaseIds[n-1]);
    
    if(gatherHaloFlags[n]!=1){
      (*nonHalo)->Ngather += test;
      ++nNonHalo;
    }
  }
  
  (*halo)->gatherOffsets  = (iint*) calloc((*halo)->Ngather+1, sizeof(iint));
  (*halo)->gatherLocalIds = (iint*) calloc(nHalo, sizeof(iint));
  (*halo)->gatherBaseIds  = (iint*) calloc((*halo)->Ngather, sizeof(iint));

  (*nonHalo)->gatherOffsets  = (iint*) calloc((*nonHalo)->Ngather+1, sizeof(iint)); 
  (*nonHalo)->gatherLocalIds = (iint*) calloc(nNonHalo, sizeof(iint));
  (*nonHalo)->gatherBaseIds  = (iint*) calloc((*nonHalo)->Ngather, sizeof(iint));
  
  // only finds bases
  nHalo = 0;
  nNonHalo = 0;
  (*halo)->Ngather = 0; // reset counter
  (*nonHalo)->Ngather = 0; // reset counter

#if 0
  for(iint n=0;n<Nlocal;++n){
    printf("rank%d: n=%d, base=%d, local=%d, halo=%d\n", rank, n, gatherBaseIds[n], gatherLocalIds[n], gatherHaloFlags[n]);
  }
#endif
  
  for(iint n=0;n<Nlocal;++n){
    iint test = (n==0) ? 1: (gatherBaseIds[n] != gatherBaseIds[n-1]);

    // increment unique base counter and record index into shuffled list of nodes
    if(gatherHaloFlags[n]==1){
      if(test){
        (*halo)->gatherOffsets[(*halo)->Ngather] = nHalo;  
        (*halo)->gatherBaseIds[(*halo)->Ngather] = gatherBaseIds[n];
  ++((*halo)->Ngather);
      }
      (*halo)->gatherLocalIds[nHalo] = gatherLocalIds[n];
      ++nHalo;
    }
  }
  
  for(iint n=0;n<Nlocal;++n){

    iint test = (n==0) ? 1: (gatherBaseIds[n] != gatherBaseIds[n-1]);
    
    if(gatherHaloFlags[n]!=1){
      if(test){
  (*nonHalo)->gatherOffsets[(*nonHalo)->Ngather] = nNonHalo;
  ++((*nonHalo)->Ngather);
      }
      (*nonHalo)->gatherLocalIds[nNonHalo] = gatherLocalIds[n];
      ++nNonHalo;
    }
  }
  (*halo)->gatherOffsets[(*halo)->Ngather] = nHalo;
  (*nonHalo)->gatherOffsets[(*nonHalo)->Ngather] = nNonHalo;
    
  // if there are halo nodes to gather
  if((*halo)->Ngather){

    occa::memory o_gatherTmpPinned = mesh->device.mappedAlloc((*halo)->Ngather*Nbytes, NULL);
    (*halo)->gatherTmp = (char*) o_gatherTmpPinned.getMappedPointer(); // (char*) calloc((*halo)->Ngather*Nbytes, sizeof(char));
    printf("host gatherTmp = %p, Nbytes = %d, Ngather = %d\n",  o_gatherTmpPinned.getMappedPointer(), Nbytes, (*halo)->Ngather);

    //    (*halo)->gatherTmp = (char*) calloc((*halo)->Ngather*Nbytes, sizeof(char));
    
    (*halo)->o_gatherTmp      = mesh->device.malloc((*halo)->Ngather*Nbytes,           (*halo)->gatherTmp);
    (*halo)->o_gatherOffsets  = mesh->device.malloc(((*halo)->Ngather+1)*sizeof(iint), (*halo)->gatherOffsets);
    (*halo)->o_gatherLocalIds = mesh->device.malloc(nHalo*sizeof(iint),                (*halo)->gatherLocalIds);
    
    // initiate gslib gather-scatter comm pattern on halo nodes only
    (*halo)->gatherGsh = gsParallelGatherScatterSetup((*halo)->Ngather, (*halo)->gatherBaseIds);
  }

  // if there are non-halo nodes to gather
  if((*nonHalo)->Ngather){

    (*nonHalo)->gatherGsh = NULL;
  
    (*nonHalo)->o_gatherOffsets  = mesh->device.malloc(((*nonHalo)->Ngather+1)*sizeof(iint), (*nonHalo)->gatherOffsets);
    (*nonHalo)->o_gatherLocalIds = mesh->device.malloc(nNonHalo*sizeof(iint),                (*nonHalo)->gatherLocalIds);
  }
  return;
}
