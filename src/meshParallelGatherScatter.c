#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh.h"

// ok to use o_v = o_gsv
void meshParallelGatherScatter(mesh_t *mesh,
                               ogs_t *ogs, 
                               occa::memory &o_v,
                               occa::memory &o_gsv){

  // gather on DEVICE
  mesh->gatherKernel(ogs->NtotalGather, ogs->o_gatherOffsets, ogs->o_gatherLocalIds, o_v, ogs->o_gatherTmp);

  // extract gathered halo node data [i.e. shared nodes ]
  if(ogs->Nhalo){
    mesh->getKernel(ogs->Nhalo, ogs->o_gatherTmp, ogs->o_haloLocalIds, ogs->o_haloTmp); // subv = v[ids]
    
    // copy partially gathered halo data from DEVICE to HOST
    ogs->o_haloTmp.copyTo(ogs->haloTmp);
    
    // gather across MPI processes then scatter back
    gsParallelGatherScatter(ogs->haloGsh, ogs->haloTmp, dfloatString, "add"); // warning: hardwired dfloat type
    
    // copy totally gather halo data back from HOST to DEVICE
    ogs->o_haloTmp.copyFrom(ogs->haloTmp);
    
    // insert totally gathered halo node data - need this kernel 
    mesh->putKernel(ogs->Nhalo, ogs->o_haloTmp,ogs->o_haloLocalIds, ogs->o_gatherTmp); // v[ids] = subv
  }

  // do scatter back to local nodes
  mesh->scatterKernel(ogs->NtotalGather, ogs->o_gatherOffsets, ogs->o_gatherLocalIds, ogs->o_gatherTmp, o_gsv);
}

void meshParallelGather(mesh_t *mesh,
                        ogs_t *ogs, 
                        occa::memory &o_v,
                        occa::memory &o_gv){

  // MPI info
  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(ogs->Nhalo){
    // gather on DEVICE
    mesh->gatherKernel(ogs->NtotalGather, ogs->o_gatherOffsets, ogs->o_gatherLocalIds, o_v, ogs->o_gatherTmp);
    
    mesh->getKernel(ogs->Nhalo, ogs->o_gatherTmp, ogs->o_haloLocalIds, ogs->o_haloTmp); 

    // copy partially gathered halo data from DEVICE to HOST
    ogs->o_haloTmp.copyTo(ogs->haloTmp);
    
    // gather across MPI processes then scatter back
    gsParallelGatherScatter(ogs->haloGsh, ogs->haloTmp, dfloatString, "add"); // warning: hardwired dfloat type adn op
    
    if (ogs->NownedHalo) {
      // copy owned and gathered halo data back from HOST to DEVICE
      ogs->o_haloTmp.copyFrom(ogs->haloTmp, ogs->NownedHalo*sizeof(dfloat));
    
      // insert totally gathered halo node data - need this kernel 
      mesh->putKernel(ogs->NownedHalo, ogs->o_haloTmp, ogs->o_haloLocalIds, ogs->o_gatherTmp); 
    }
    o_gv.copyFrom(ogs->o_gatherTmp,ogs->Ngather*sizeof(dfloat)); //copy the gathered data
  } else {
    // gather on DEVICE
    mesh->gatherKernel(ogs->Ngather, ogs->o_gatherOffsets, ogs->o_gatherLocalIds, o_v, o_gv);
  }
}

void meshParallelScatter(mesh_t *mesh,
                        ogs_t *ogs, 
                        occa::memory &o_v,
                        occa::memory &o_sv){

  // MPI info
  iint rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(ogs->Nhalo){
    if (ogs->NownedHalo)
      mesh->getKernel(ogs->NownedHalo, o_v, ogs->o_haloLocalIds, ogs->o_haloTmp); 
    
    //zero out excess halo buffer
    for (iint n=ogs->NownedHalo;n<ogs->Nhalo;n++) ogs->haloTmp[n] = 0;

    if (ogs->NownedHalo)  // copy owned halo data from DEVICE to HOST
      ogs->o_haloTmp.copyTo(ogs->haloTmp, ogs->NownedHalo*sizeof(dfloat));
    
    o_v.copyTo(ogs->o_gatherTmp,ogs->Ngather*sizeof(dfloat)); //copy the data

    // scatter across MPI processes 
    gsParallelGatherScatter(ogs->haloGsh, ogs->haloTmp, dfloatString, "add"); // warning: hardwired dfloat type adn op
    
    // copy scattered halo data from DEVICE to HOST
    ogs->o_haloTmp.copyFrom(ogs->haloTmp);

    mesh->putKernel(ogs->Nhalo, ogs->o_haloTmp, ogs->o_haloLocalIds, ogs->o_gatherTmp);

    // do scatter back to local nodes
    mesh->scatterKernel(ogs->NtotalGather, ogs->o_gatherOffsets, ogs->o_gatherLocalIds, ogs->o_gatherTmp, o_sv);
  } else {
    // do scatter back to local nodes
    mesh->scatterKernel(ogs->Ngather, ogs->o_gatherOffsets, ogs->o_gatherLocalIds, o_v, o_sv);  
  }
}