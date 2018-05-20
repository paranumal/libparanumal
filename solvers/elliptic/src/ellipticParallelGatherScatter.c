#include "elliptic.h"

void ellipticParallelGatherScatter(mesh_t *mesh, ogs_t *ogs, occa::memory &o_q, const char *type, const char *op){

  // use gather map for gather and scatter
  occaTimerTic(mesh->device,"meshParallelGatherScatter");
  meshParallelGatherScatter(mesh, ogs, o_q);
  occaTimerToc(mesh->device,"meshParallelGatherScatter");
}
