#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh3D.h"

typedef struct {

  occa::memory o_vmapPP;
  occa::memory o_faceNodesP;
  occa::memory o_oasForward;
  occa::memory o_oasBack;
  occa::memory o_oasDiagInvOp;

  occa::kernel restrictKernel;
  occa::kernel preconKernel;
  
  ogs_t *ogsP;

  occa::memory o_diagA;
  
} precon_t;

void ellipticRunHex3D(mesh3D *mesh);

void ellipticOccaRunHex3D(mesh3D *mesh);

void ellipticSetupHex3D(mesh3D *mesh, ogs_t **ogs, precon_t **precon, dfloat lambda);

void ellipticVolumeHex3D(mesh3D *mesh);

void ellipticSurfaceHex3D(mesh3D *mesh, dfloat time);

void ellipticUpdateHex3D(mesh3D *mesh, dfloat rka, dfloat rkb);

void ellipticErrorHex3D(mesh3D *mesh, dfloat time);

void ellipticParallelGatherScatter3D(mesh3D *mesh, ogs_t *ogs, occa::memory &o_v, occa::memory &o_gsv, const char *type);

precon_t *ellipticOasPreconSetupHex3D(mesh3D *mesh, ogs_t *ogs, dfloat lambda);

