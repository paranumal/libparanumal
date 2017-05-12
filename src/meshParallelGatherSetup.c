#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include <mpi.h>

#include "mesh.h"

typedef struct {

  iint localId;
  iint globalId;
  iint ownerRank;
  
}parallelNode_t;

// assume nodes locally sorted by rank then global index
// assume gather and scatter are the same sets
hgs_t *meshParallelGatherSetup(mesh_t *mesh,    // provides DEVICE
				      iint Nlocal,     // number of local nodes
				      iint Nbytes,     // number of bytes per node
				      iint *globalNumbering,  // global index of nodes
				      iint *globalOwners){   // node owners

  iint rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  hgs_t *hgs = (hgs_t*) calloc(1, sizeof(hgs_t));

  hgs->Nlocal = Nlocal;

  parallelNode_t *nodes = (parallelNode_t *) calloc(Nlocal,sizeof(parallelNode_t));

  for (iint n=0;n<Nlocal;n++) {
    nodes[n].localId = n;
    nodes[n].globalId = globalNumbering[n];
    nodes[n].ownerRank = globalOwners[n];
  }

  //sort by owner rank
  qsort(nodes, Nlocal, sizeof(parallelNode_t), parallelCompareOwners);

  hgs->Nsend = (iint *) calloc(size,sizeof(iint));
  for (iint n=0;n<Nlocal;n++)    
    hgs->Nsend[nodes[n].ownerRank]++

  //move the nodes that are already on the right rank to the front
  iint cnt =0;
  for (iint n=0;n<Nlocal;n++) {
    if (nodes[n].ownerRank == rank) {
      parallelNode_t tmp = nodes[cnt];
      nodes[cnt++] = nodes[n];
      nodes[n] = tmp;
    }
  }

  //record this ordering
  hgs->sortIds = (iint *) calloc(Nlocal,sizeof(iint));
  for (iint n=0;n<Nlocal;n++)
    hgs->sortIds[nodes[n].localId] = n;

  //share how many nodes we're sending
  hgs->Nrecv = (iint *) calloc(size,sizeof(iint));
  MPI_Alltoall(hgs->Nsend, 1, MPI_IINT, hgs->Nrecv, 1, MPI_IINT, MPI_COMM_WORLD);

  

  // ------------------------------------------------------------  
  // 0. propagate halo flags uniformly using a disposable gs instance



  void *allGsh = gsParallelGatherScatterSetup(Nlocal, gatherBaseIds);

  // compute max halo flag using global numbering
  gsParallelGatherScatter(allGsh, gatherHaloFlags, "int", "max"); // should use iint

  // tidy up
  gsParallelGatherScatterDestroy(allGsh);

  // ------------------------------------------------------------
  // 1. count number of unique base nodes on this process
  ogs->Ngather = 0; // assumes at least one base node
  for(iint n=0;n<Nlocal;++n){
    iint test = (n==0) ? 1: (gatherBaseIds[n] != gatherBaseIds[n-1]);
    ogs->Ngather += test;
  }
  
  iint *gatherOffsets = (iint*) calloc(ogs->Ngather+1, sizeof(iint)); // offset into sorted list of nodes

  // only finds bases
  ogs->Ngather = 0; // reset counter
  for(iint n=0;n<Nlocal;++n){
    iint test = (n==0) ? 1: (gatherBaseIds[n] != gatherBaseIds[n-1]);
    if(test){
      gatherOffsets[ogs->Ngather++] = n;  // increment unique base counter and record index into shuffled list of ndoes
    }
  }
  gatherOffsets[ogs->Ngather] = Nlocal;

  char *gatherTmp = (char*) calloc(ogs->Ngather*Nbytes, sizeof(char));
  
  // allocate buffers on DEVICE
  ogs->o_gatherTmp      = mesh->device.malloc(ogs->Ngather*Nbytes, gatherTmp);
  ogs->o_gatherOffsets  = mesh->device.malloc((ogs->Ngather+1)*sizeof(iint), gatherOffsets);
  ogs->o_gatherLocalIds = mesh->device.malloc(Nlocal*sizeof(iint), gatherLocalIds);
  
  ogs->gatherLocalIds  = gatherLocalIds;
  ogs->gatherBaseIds   = gatherBaseIds;
  ogs->gatherHaloFlags = gatherHaloFlags;

  // list of nodes to extract from DEVICE gathered array
  ogs->Nhalo = 0; 
  for(iint n=0;n<ogs->Ngather;++n){ // could be this? 
    for(iint id=gatherOffsets[n];id<gatherOffsets[n+1];++id){
      if(gatherHaloFlags[id]){ // if any of these are labelled as halo then mark
	++ogs->Nhalo;
	break;
      }
    }
  }
  
  // set up gather-scatter of halo nodes
  ogs->gatherGsh = NULL;
  if(ogs->Nhalo){
    
    iint *haloLocalIds  = (iint*) calloc(ogs->Nhalo, sizeof(iint));
    iint *haloGlobalIds = (iint*) calloc(ogs->Nhalo, sizeof(iint));

    ogs->Nhalo = 0;
    for(iint n=0;n<ogs->Ngather;++n){
      for(iint id=gatherOffsets[n];id<gatherOffsets[n+1];++id){
      	if(gatherHaloFlags[id]){
      	  haloLocalIds[ogs->Nhalo] = n;
      	  haloGlobalIds[ogs->Nhalo] = gatherBaseIds[id];
      	  ++ogs->Nhalo;
      	  break;
      	}
      }
    }     
    
    // list of halo node indices
    ogs->o_haloLocalIds = mesh->device.malloc(ogs->Nhalo*sizeof(iint), haloLocalIds);
    
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
