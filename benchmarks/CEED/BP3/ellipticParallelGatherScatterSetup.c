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
  
  // ------------------------------------------------------------  
  // 0. propagate halo flags uniformly using a disposable gs instance
  void *allGsh = gsParallelGatherScatterSetup(Nlocal, gatherBaseIds);

  // compute max halo flag using global numbering
  gsParallelGatherScatter(allGsh, gatherHaloFlags, "int", "max"); // should use iint

  // tidy up
  gsParallelGatherScatterDestroy(allGsh);

  // initialize gather structs
  *halo = (ogs_t*) calloc(1, sizeof(ogs_t));
  *nonHalo = (ogs_t*) calloc(1, sizeof(ogs_t));
  
  // ------------------------------------------------------------
  // 1. count number of unique base nodes on this process 
  (*halo)->Ngather = 0;
  (*nonHalo)->Ngather = 0;

  iint NhaloLocal = 0;
  iint NnonHaloLocal = 0;
  
  for(iint n=0;n<Nlocal;++n){
    iint test = (n==0) ? 1: (gatherBaseIds[n] != gatherBaseIds[n-1]);
    if(gatherHaloFlags[n]==1)
      (*halo)->Ngather += test;
    else
      (*nonHalo)->Ngather += test;
    
    // total number of halo nodes
    NhaloLocal += (gatherHaloFlags[n]==1);
    NnonHaloLocal += (gatherHaloFlags[n]!=1);
  }
  
  (*halo)->gatherOffsets    = (iint*) calloc((*halo)->Ngather+1, sizeof(iint)); 
  (*nonHalo)->gatherOffsets = (iint*) calloc((*nonHalo)->Ngather+1, sizeof(iint)); 

  (*halo)->gatherLocalIds = (iint*) calloc(NhaloLocal, sizeof(iint));
  (*nonHalo)->gatherLocalIds = (iint*) calloc(NnonHaloLocal, sizeof(iint));

  (*halo)->gatherBaseIds = (iint*) calloc(NhaloLocal, sizeof(iint));
  (*nonHalo)->gatherBaseIds = (iint*) calloc(NnonHaloLocal, sizeof(iint));
  
  // only finds bases
  iint nHalo = 0;
  iint nNonHalo = 0;
  (*halo)->Ngather = 0; // reset counter
  (*nonHalo)->Ngather = 0; // reset counter

  for(iint n=0;n<Nlocal;++n){
    iint test = (n==0) ? 1: (gatherBaseIds[n] != gatherBaseIds[n-1]);
    if(test){
      // increment unique base counter and record index into shuffled list of ndoes
      if(gatherHaloFlags[n]==1)
	(*halo)->gatherOffsets[(*halo)->Ngather++] = nHalo;  
      else
	(*nonHalo)->gatherOffsets[(*nonHalo)->Ngather++] = nNonHalo;  
    }
    
    if(gatherHaloFlags[n]==1){
      (*halo)->gatherLocalIds[nHalo] = n;
      (*halo)->gatherBaseIds[nHalo++] = gatherBaseIds[n];
    }else{
      (*nonHalo)->gatherLocalIds[nNonHalo] = n;
      (*nonHalo)->gatherBaseIds[nNonHalo++] = gatherBaseIds[n];
    }
  }
  
  (*halo)->gatherOffsets[(*halo)->Ngather]          = (*halo)->Ngather + 1;
  (*nonHalo)->gatherOffsets[(*nonHalo)->Ngather] = (*nonHalo)->Ngather + 1;

  (*halo)->gatherTmp = (char*) calloc((*halo)->Ngather*Nbytes, sizeof(char));
  (*nonHalo)->gatherTmp = (char*) calloc((*nonHalo)->Ngather*Nbytes, sizeof(char));

  (*halo)->gatherGsh = NULL;  
  (*nonHalo)->gatherGsh = NULL;
  
  // allocate buffers on DEVICE
  if((*halo)->Ngather){
    (*halo)->o_gatherTmp      = mesh->device.malloc((*halo)->Ngather*Nbytes,           (*halo)->gatherTmp);
    (*halo)->o_gatherOffsets  = mesh->device.malloc(((*halo)->Ngather+1)*sizeof(iint), (*halo)->gatherOffsets);
    (*halo)->o_gatherLocalIds = mesh->device.malloc((*halo)->Ngather*sizeof(iint),     (*halo)->gatherLocalIds);

    // initiate gslib gather-scatter comm pattern on halo nodes only
    (*halo)->gatherGsh = gsParallelGatherScatterSetup((*halo)->Ngather, (*halo)->gatherBaseIds);
  }

  if((*nonHalo)->Ngather){
    (*nonHalo)->o_gatherTmp      = mesh->device.malloc((*nonHalo)->Ngather*Nbytes,           (*nonHalo)->gatherTmp);
    (*nonHalo)->o_gatherOffsets  = mesh->device.malloc(((*nonHalo)->Ngather+1)*sizeof(iint), (*nonHalo)->gatherOffsets);
    (*nonHalo)->o_gatherLocalIds = mesh->device.malloc((*nonHalo)->Ngather*sizeof(iint),     (*nonHalo)->gatherLocalIds);
  }
  return;
}
