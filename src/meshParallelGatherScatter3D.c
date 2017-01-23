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
void meshParallelGatherScatter3D(mesh3D *mesh,
				 iint Ngather,                 // input: number of gather nodes 
				 occa::memory &o_gatherOffsets, // input: start of local bases
				 occa::memory &o_gatherLocalIds,// input: base connected nodes
				 occa::memory &o_gatherTmp,     // input: DEVICE gather buffer
				 iint Nhalo,                   // input: number of halo nodes
				 occa::memory &o_haloLocalIds,  // input: list of halo nodes to
				 occa::memory &o_haloTmp,       // input: temporary halo buffer
				 void *haloTmp,                // input: temporary HOST halo buffer
				 iint Nscatter,                 // input: number of scatter nodes 
				 occa::memory &o_scatterOffsets, // input: start of local bases
				 occa::memory &o_scatterLocalIds,// input: base connected nodes
				 void *gsh,
				 occa::memory &o_v,
				 occa::memory &o_gsv,
				 const char *type){

  // gather on DEVICE
  mesh->gatherKernel(Ngather, o_gatherOffsets, o_gatherLocalIds, o_v, o_gatherTmp);

  // extract gathered halo node data [i.e. shared nodes ]
  if(Nhalo){
    mesh->getKernel(Nhalo, o_gatherTmp, o_haloLocalIds, o_haloTmp); // subv = v[ids]
    
    // copy partially gathered halo data from DEVICE to HOST
    o_haloTmp.copyTo(haloTmp);
    
    // gather across MPI processes then scatter back
    gsParallelGatherScatter(gsh, haloTmp, dfloatString); // danger on hardwired type
    
    // copy totally gather halo data back from HOST to DEVICE
    o_haloTmp.copyFrom(haloTmp);
    
    // insert totally gathered halo node data - need this kernel 
    mesh->putKernel(Nhalo, o_haloTmp,o_haloLocalIds, o_gatherTmp); // v[ids] = subv
  }
  
  // do scatter back to local nodes
  mesh->scatterKernel(Nscatter, o_scatterOffsets, o_scatterLocalIds, o_gatherTmp, o_gsv);
  
}
