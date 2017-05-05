#include "ellipticHex3D.h"

void ellipticParallelGatherScatterHex3D(mesh3D *mesh, ogs_t *ogs, occa::memory &o_q, occa::memory &o_gsq, const char *type, const char *op){

  mesh->device.finish();
  occa::tic("meshParallelGatherScatter3D");
  
  // use gather map for gather and scatter
  meshParallelGatherScatter(mesh, ogs, o_q, o_gsq, type, op);

  mesh->device.finish();
  occa::toc("meshParallelGatherScatter3D");
  
}