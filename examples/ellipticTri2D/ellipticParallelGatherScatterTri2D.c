#include "ellipticTri2D.h"

void ellipticParallelGatherScatterTri2D(mesh2D *mesh, ogs_t *ogs, occa::memory &o_q, occa::memory &o_gsq, const char *type, const char *op){

  // use gather map for gather and scatter
  occaTimerTic(mesh->device,"meshParallelGatherScatter2D");
  meshParallelGatherScatter(mesh, ogs, o_q, o_gsq);
  occaTimerToc(mesh->device,"meshParallelGatherScatter2D");
}
