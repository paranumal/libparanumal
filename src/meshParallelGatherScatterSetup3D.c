#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include <mpi.h>

#include "mesh3D.h"

extern "C"
{
  void *gsParallelGatherScatterSetup(int Ngather, int *gatherIds);
}

// assume nodes locally sorted by rank then global index
void meshParallelGatherScatterSetup3D(mesh3D *mesh,    // provides DEVICE
				      iint Nlocal,     // number of local nodes
				      iint Nbytes,     // number of bytes per node
				      iint *localIds,  // local index of nodes
				      iint *baseIds,   // global index of their base nodes
				      iint *baseRanks, // rank of their base nodes
				      iint *maxRanks,  // max rank connected to base node
				      iint &Ngather,   // output: number of gather nodes 
				      occa::memory &o_gatherOffsets, // output: start of local bases
				      occa::memory &o_gatherLocalIds,// output: base connected nodes
				      occa::memory &o_gatherTmp,     // output: DEVICE gather buffer
				      iint &Nhalo,     // output: number of halo nodes
				      occa::memory &o_haloLocalIds, // list of halo nodes to
				      occa::memory &o_haloTmp, // temporary halo buffer
				      void **haloTmp, // temporary HOST halo buffer
				      void **gsh){ // output: gather-scatter
  
  // 1. count number of unique base nodes on this process
  Ngather = 0; // assumes at least one base node
  for(iint n=0;n<Nlocal;++n){
    iint test = (n==0) ? 1: (baseIds[n] != baseIds[n-1]);
    Ngather += test;
  }
  
  iint *gatherIds     = (iint*) calloc(Ngather, sizeof(iint)); // global labels for unique nodes
  iint *gatherOffsets = (iint*) calloc(Ngather+1, sizeof(iint)); // offset into sorted list of nodes

  // only finds bases
  Ngather = 0; // reset counter
  for(iint n=0;n<Nlocal;++n){
    iint test = (n==0) ? 1: (baseIds[n] != baseIds[n-1]);
    if(test){
      gatherIds[Ngather] = baseIds[n]; // global indices of base nodes
      gatherOffsets[Ngather++] = n;  // increment unique base counter and record index into shuffled list of ndoes
    }
  }
  gatherOffsets[Ngather] = Nlocal;

  // allocate buffers on DEVICE
  o_gatherTmp      = mesh->device.malloc(Ngather*Nbytes);
  o_gatherOffsets  = mesh->device.malloc((Ngather+1)*sizeof(iint), gatherOffsets);
  o_gatherLocalIds = mesh->device.malloc(Nlocal*sizeof(iint), localIds);
  
  // list of nodes to extract from DEVICE gathered array
  Nhalo = 0;
  for(iint n=0;n<Ngather;++n){
    int id = gatherOffsets[n];
    if(baseRanks[id]!=maxRanks[id]){
      ++Nhalo;
    }
  }
  
  printf("Nhalo = %d, Ngather = %d\n", Nhalo, Ngather);

  // set up gather-scatter of halo nodes
  *gsh = NULL;
  if(Nhalo){
    
    iint *haloLocalIds = (iint*) calloc(Nhalo, sizeof(iint));
    iint *haloGlobalIds = (iint*) calloc(Nhalo, sizeof(iint));
    Nhalo = 0;
    for(iint n=0;n<Ngather;++n){
      int id = gatherOffsets[n];
      if(baseRanks[id]!=maxRanks[id]){
	haloLocalIds[Nhalo] = n;
	haloGlobalIds[Nhalo] = baseIds[id];
	++Nhalo;
      }
    }     
    
    // list of halo node indices
    o_haloLocalIds = mesh->device.malloc(Nhalo*sizeof(iint), haloLocalIds);
    
    // allocate buffer for gathering on halo nodes (danger on size of buffer)
    *haloTmp = (dfloat*) malloc(Nhalo*Nbytes);
    o_haloTmp  = mesh->device.malloc(Nhalo*Nbytes, *haloTmp);
    
    // initiate gslib gather-scatter comm pattern
    *gsh = gsParallelGatherScatterSetup(Nhalo, haloGlobalIds);
    
    free(haloGlobalIds);
    free(haloLocalIds);
  }
  
  free(gatherIds);
  free(gatherOffsets);
}
