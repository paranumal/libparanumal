#include "ellipticTri2D.h"

void ellipticParallelGatherScatterTri2D(mesh2D *mesh, ogs_t *ogs, occa::memory &o_q, occa::memory &o_gsq, const char *type, const char *op){

  mesh->device.finish();
  occa::tic("meshParallelGatherScatter2D");
  
  // use gather map for gather and scatter
  meshParallelGatherScatter(mesh, ogs, o_q, o_gsq, type, op);

  mesh->device.finish();
  occa::toc("meshParallelGatherScatter2D");
  
}
