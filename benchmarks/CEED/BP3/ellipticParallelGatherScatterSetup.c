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

#include "ellipticHex3D.h"

// assume nodes locally sorted by rank then global index
// assume gather and scatter are the same sets
void ellipticParallelGatherScatterSetup(mesh_t *mesh,    // provides DEVICE
                                        int Nlocal,     // number of local nodes
                                        int Nbytes,     // number of bytes per node
                                        int *gatherLocalIds,  // local index of nodes
                                        int *gatherBaseIds,   // global index of their base nodes
                                        int *gatherHaloFlags,
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
  gsParallelGatherScatter(allGsh, gatherHaloFlags, "int", "max"); // should use int

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

  int nHalo = 0;
  int nNonHalo = 0;
  
  for(int n=0;n<Nlocal;++n){
    int test = (n==0) ? 1: (gatherBaseIds[n] != gatherBaseIds[n-1]);
    if(gatherHaloFlags[n]==1){
      (*halo)->Ngather += test;
      ++nHalo;
    }
  }


  for(int n=0;n<Nlocal;++n){
    int test = (n==0) ? 1: (gatherBaseIds[n] != gatherBaseIds[n-1]);
    
    if(gatherHaloFlags[n]!=1){
      (*nonHalo)->Ngather += test;
      ++nNonHalo;
    }
  }
  
  (*halo)->gatherOffsets  = (int*) calloc((*halo)->Ngather+1, sizeof(int));
  (*halo)->gatherLocalIds = (int*) calloc(nHalo, sizeof(int));
  (*halo)->gatherBaseIds  = (int*) calloc((*halo)->Ngather, sizeof(int));

  (*nonHalo)->gatherOffsets  = (int*) calloc((*nonHalo)->Ngather+1, sizeof(int)); 
  (*nonHalo)->gatherLocalIds = (int*) calloc(nNonHalo, sizeof(int));
  (*nonHalo)->gatherBaseIds  = (int*) calloc((*nonHalo)->Ngather, sizeof(int));
  
  // only finds bases
  nHalo = 0;
  nNonHalo = 0;
  (*halo)->Ngather = 0; // reset counter
  (*nonHalo)->Ngather = 0; // reset counter

#if 0
  for(int n=0;n<Nlocal;++n){
    printf("rank%d: n=%d, base=%d, local=%d, halo=%d\n", rank, n, gatherBaseIds[n], gatherLocalIds[n], gatherHaloFlags[n]);
  }
#endif
  
  for(int n=0;n<Nlocal;++n){
    int test = (n==0) ? 1: (gatherBaseIds[n] != gatherBaseIds[n-1]);

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
  
  for(int n=0;n<Nlocal;++n){

    int test = (n==0) ? 1: (gatherBaseIds[n] != gatherBaseIds[n-1]);
    
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
    //    printf("host gatherTmp = %p, Nbytes = %d, Ngather = %d\n",  o_gatherTmpPinned.getMappedPointer(), Nbytes, (*halo)->Ngather);

    //    (*halo)->gatherTmp = (char*) calloc((*halo)->Ngather*Nbytes, sizeof(char));
    
    (*halo)->o_gatherTmp      = mesh->device.malloc((*halo)->Ngather*Nbytes,           (*halo)->gatherTmp);
    (*halo)->o_gatherOffsets  = mesh->device.malloc(((*halo)->Ngather+1)*sizeof(int), (*halo)->gatherOffsets);
    (*halo)->o_gatherLocalIds = mesh->device.malloc(nHalo*sizeof(int),                (*halo)->gatherLocalIds);
    
    // initiate gslib gather-scatter comm pattern on halo nodes only
    (*halo)->gatherGsh = gsParallelGatherScatterSetup((*halo)->Ngather, (*halo)->gatherBaseIds);
  }

  // if there are non-halo nodes to gather
  if((*nonHalo)->Ngather){

    (*nonHalo)->gatherGsh = NULL;
  
    (*nonHalo)->o_gatherOffsets  = mesh->device.malloc(((*nonHalo)->Ngather+1)*sizeof(int), (*nonHalo)->gatherOffsets);
    (*nonHalo)->o_gatherLocalIds = mesh->device.malloc(nNonHalo*sizeof(int),                (*nonHalo)->gatherLocalIds);
  }
  return;
}
