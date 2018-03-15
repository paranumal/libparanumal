#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh.h"

// assume nodes locally sorted by rank then global index
// assume gather and scatter are the same sets
ogs_t *meshParallelGatherScatterSetup(mesh_t *mesh,    // provides DEVICE
				      int Nlocal,     // number of local nodes
				      int Nbytes,     // number of bytes per node
				      int *gatherLocalIds,  // local index of nodes
				      int *gatherBaseIds,   // global index of their base nodes
				      int *gatherHaloFlags){   // 1 for halo node, 0 for not

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  ogs_t *ogs = (ogs_t*) calloc(1, sizeof(ogs_t));

  // ------------------------------------------------------------
  // 0. propagate halo flags uniformly using a disposable gs instance

  void *allGsh = gsParallelGatherScatterSetup(Nlocal, gatherBaseIds);

  // compute max halo flag using global numbering
  gsParallelGatherScatter(allGsh, gatherHaloFlags, "int", "max"); // should use int

  // tidy up
  gsParallelGatherScatterDestroy(allGsh);

  // ------------------------------------------------------------
  // 1. count number of unique base nodes on this process
  ogs->Ngather = 0; // assumes at least one base node
  for(int n=0;n<Nlocal;++n){
    int test = (n==0) ? 1: (gatherBaseIds[n] != gatherBaseIds[n-1]);
    ogs->Ngather += test;
  }

  int *gatherOffsets = (int*) calloc(ogs->Ngather+1, sizeof(int)); // offset into sorted list of nodes

  // only finds bases
  ogs->Ngather = 0; // reset counter
  for(int n=0;n<Nlocal;++n){
    int test = (n==0) ? 1: (gatherBaseIds[n] != gatherBaseIds[n-1]);
    if(test){
      gatherOffsets[ogs->Ngather++] = n;  // increment unique base counter and record index into shuffled list of ndoes
    }
  }
  gatherOffsets[ogs->Ngather] = Nlocal;

  char *gatherTmp = (char*) calloc(ogs->Ngather*Nbytes, sizeof(char));

  // allocate buffers on DEVICE
  ogs->o_gatherTmp      = mesh->device.malloc(ogs->Ngather*Nbytes, gatherTmp);
  ogs->o_gatherOffsets  = mesh->device.malloc((ogs->Ngather+1)*sizeof(int), gatherOffsets);
  ogs->o_gatherLocalIds = mesh->device.malloc(Nlocal*sizeof(int), gatherLocalIds);

  ogs->gatherLocalIds  = gatherLocalIds;
  ogs->gatherBaseIds   = gatherBaseIds;
  ogs->gatherHaloFlags = gatherHaloFlags;

  // list of nodes to extract from DEVICE gathered array
  ogs->Nhalo = 0;
  for(int n=0;n<ogs->Ngather;++n){ // could be this?
    for(int id=gatherOffsets[n];id<gatherOffsets[n+1];++id){
      if(gatherHaloFlags[id]){ // if any of these are labelled as halo then mark
	++ogs->Nhalo;
	break;
      }
    }
  }

  // set up gather-scatter of halo nodes
  ogs->gatherGsh = NULL;
  if(ogs->Nhalo){

    int *haloLocalIds  = (int*) calloc(ogs->Nhalo, sizeof(int));
    int *haloGlobalIds = (int*) calloc(ogs->Nhalo, sizeof(int));

    ogs->Nhalo = 0;
    for(int n=0;n<ogs->Ngather;++n){
      for(int id=gatherOffsets[n];id<gatherOffsets[n+1];++id){
      	if(gatherHaloFlags[id]){
      	  haloLocalIds[ogs->Nhalo] = n;
      	  haloGlobalIds[ogs->Nhalo] = gatherBaseIds[id];
      	  ++ogs->Nhalo;
      	  break;
      	}
      }
    }

    // list of halo node indices
    ogs->o_haloLocalIds = mesh->device.malloc(ogs->Nhalo*sizeof(int), haloLocalIds);

    // allocate buffer for gathering on halo nodes (danger on size of buffer)
    if(Nbytes != sizeof(dfloat)) printf("DANGER WILL ROBINSON\n");

    ogs->haloTmp = (dfloat*) calloc(ogs->Nhalo, sizeof(dfloat));
    ogs->o_haloTmp  = mesh->device.malloc(ogs->Nhalo*Nbytes, ogs->haloTmp);

    // initiate gslib gather-scatter comm pattern
    ogs->gatherGsh = gsParallelGatherScatterSetup(ogs->Nhalo, haloGlobalIds);

    free(haloGlobalIds);
    free(haloLocalIds);
  }

  free(gatherOffsets);

  return ogs;
}
