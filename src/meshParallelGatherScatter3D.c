#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh3D.h"
#include "mpi.h"

extern "C"
{
  void gsParallelGatherScatter(void *gsh, void *v, const char *type);
}

// ok to use o_v = o_gsv
void meshParallelGatherScatter3D(mesh3D *mesh, occa::memory &o_v, occa::memory &o_gsv, const char *type){

  // gather on DEVICE
  mesh->gatherKernel(mesh->NuniqueBases, mesh->o_gatherNodeOffsets, mesh->o_gatherLocalNodes, o_v, mesh->o_gatherTmp);

  // extract gathered halo node data [i.e. shared nodes ]
  if(mesh->NnodeHalo){
    mesh->getKernel(mesh->NnodeHalo, mesh->o_gatherTmp, mesh->o_nodeHaloIds, mesh->o_subGatherTmp); // subv = v[ids]
    
    // copy partially gathered halo data from DEVICE to HOST
    mesh->o_subGatherTmp.copyTo(mesh->subGatherTmp);
    
    // gather across MPI processes then scatter back
    gsParallelGatherScatter(mesh->gsh, mesh->subGatherTmp, "float");
    
    // copy totally gather halo data back from HOST to DEVICE
    mesh->o_subGatherTmp.copyFrom(mesh->subGatherTmp);
    
    // insert totally gathered halo node data - need this kernel 
    mesh->putKernel(mesh->NnodeHalo, mesh->o_subGatherTmp, mesh->o_nodeHaloIds, mesh->o_gatherTmp); // v[ids] = subv
  }
  
  // do scatter back to local nodes
  mesh->scatterKernel(mesh->NuniqueBases, mesh->o_gatherNodeOffsets, mesh->o_gatherLocalNodes, mesh->o_gatherTmp, o_gsv);
  
}
