#include "mesh.h"

void meshApplyElementMatrix(mesh_t *mesh, dfloat *A, dfloat *q, dfloat *Aq) {

  dfloat *Aqn = (dfloat*) calloc(mesh->Np,sizeof(dfloat));
  for (dlong e=0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Np;n++) {
      Aqn[n] = 0;
      for (int k=0;k<mesh->Np;k++) {
        Aqn[n] += A[k+n*mesh->Np]*q[k+e*mesh->Np];
      }
    }
    for (int n=0;n<mesh->Np;n++) Aq[n+e*mesh->Np] = Aqn[n];
  }
  free(Aqn);
}