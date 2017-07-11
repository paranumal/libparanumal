#include "ellipticHex3D.h"

void ellipticParallelGatherScatter(mesh3D *mesh, ogs_t *ogs, occa::memory &o_q, occa::memory &o_gsq, const char *type, const char *op){

  // use gather map for gather and scatter
  occaTimerTic(mesh->device,"meshParallelGatherScatter3D");
  meshParallelGatherScatter(mesh, ogs, o_q, o_gsq, type, op);
  occaTimerToc(mesh->device,"meshParallelGatherScatter3D");
}


void ellipticHaloGatherScatter(solver_t *solver, 
			       ogs_t *halo, 
			       occa::memory &o_v,
			       const char *type,
			       const char *op){


  mesh3D *mesh = solver->mesh;

  if(halo->Ngather){
    // set stream to halo stream
    //    mesh->device.setStream(solver->dataStream);
    
    // gather halo nodes on DEVICE
    mesh->gatherKernel(halo->Ngather, halo->o_gatherOffsets, halo->o_gatherLocalIds, o_v, halo->o_gatherTmp);
    
    // copy partially gathered halo data from DEVICE to HOST
    halo->o_gatherTmp.copyTo(halo->gatherTmp); // do I need async here (if so - need to pin gatherTmp) ?
    
    // gather across MPI processes then scatter back
    gsParallelGatherScatter(halo->gatherGsh, halo->gatherTmp, dfloatString, op); // danger on hardwired type
    
    // copy totally gather halo data back from HOST to DEVICE
    halo->o_gatherTmp.copyFrom(halo->gatherTmp); // do I need async here ?
    
    // do scatter back to local nodes
    mesh->scatterKernel(halo->Ngather, halo->o_gatherOffsets, halo->o_gatherLocalIds, halo->o_gatherTmp, o_v);

    // revert back to default stream
    //    mesh->device.setStream(solver->defaultStream);
  }
  
}


void ellipticNonHaloGatherScatter(solver_t *solver, 
				  ogs_t *nonHalo, 
				  occa::memory &o_v,
				  const char *type,
				  const char *op){


  mesh3D *mesh = solver->mesh;

  // set stream to default stream
  mesh->device.setStream(solver->defaultStream);
  
  // unified gather-scatter operation on DEVICE for non-halo nodes
  mesh->gatherScatterKernel(nonHalo->Ngather, nonHalo->o_gatherOffsets, nonHalo->o_gatherLocalIds, o_v);
  
}
