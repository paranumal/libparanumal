#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh.h"
#include "mpi.h"

extern "C"
{
  void gsParallelGatherScatter(void *gsh, void *v, const char *type, const char *op);
}

// ok to use o_v = o_gsv
void meshParallelGatherScatter(mesh_t *mesh,
			       ogs_t *ogs, 
			       occa::memory &o_v,
			       occa::memory &o_gsv,
			       const char *type,
			       const char *op){

  mesh->device.finish();
  occa::tic("gatherKernel");
  
  // gather on DEVICE
  mesh->gatherKernel(ogs->Ngather, ogs->o_gatherOffsets, ogs->o_gatherLocalIds, o_v, ogs->o_gatherTmp);

  mesh->device.finish();
  occa::toc("gatherKernel");

  //  if(strcmp(type, "float")) printf("BU55ER ; type = %s\n", type);
  
  // extract gathered halo node data [i.e. shared nodes ]
  if(ogs->Nhalo){
    mesh->getKernel(ogs->Nhalo, ogs->o_gatherTmp, ogs->o_haloLocalIds, ogs->o_haloTmp); // subv = v[ids]
    
    // copy partially gathered halo data from DEVICE to HOST
    ogs->o_haloTmp.copyTo(ogs->haloTmp);
    
    // gather across MPI processes then scatter back
    gsParallelGatherScatter(ogs->gatherGsh, ogs->haloTmp, dfloatString, op); // danger on hardwired type
    
    // copy totally gather halo data back from HOST to DEVICE
    ogs->o_haloTmp.copyFrom(ogs->haloTmp);
    
    // insert totally gathered halo node data - need this kernel 
    mesh->putKernel(ogs->Nhalo, ogs->o_haloTmp,ogs->o_haloLocalIds, ogs->o_gatherTmp); // v[ids] = subv
  }

  mesh->device.finish();
  occa::tic("scatterKernel");
  
  // do scatter back to local nodes
  //  mesh->scatterKernel(ogs->Nscatter, ogs->o_scatterOffsets, ogs->o_scatterLocalIds, ogs->o_gatherTmp, o_gsv);
  mesh->scatterKernel(ogs->Ngather, ogs->o_gatherOffsets, ogs->o_gatherLocalIds, ogs->o_gatherTmp, o_gsv);
  
  mesh->device.finish();
  occa::toc("scatterKernel");
  
}
