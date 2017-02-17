#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh2D.h"

typedef struct {

  occa::memory o_vmapPP;
  occa::memory o_faceNodesP;

  occa::memory o_oasForward;
  occa::memory o_oasBack;
  occa::memory o_oasDiagInvOp;

  occa::memory o_oasForwardDg;
  occa::memory o_oasBackDg;
  occa::memory o_oasDiagInvOpDg;

  occa::kernel restrictKernel;
  occa::kernel preconKernel;

  occa::kernel coarsenKernel;
  occa::kernel prolongateKernel;  

  ogs_t *ogsP, *ogsDg;

  occa::memory o_diagA;

  // coarse grid basis for preconditioning
  occa::memory o_V1, o_Vr1, o_Vs1, o_Vt1;
  occa::memory o_r1, o_z1;
  dfloat *r1, *z1;
  void *xxt;

  occa::memory o_coarseInvDegree;
  occa::memory o_ztmp;
  
} precon_t;

void ellipticRunTri2D(mesh2D *mesh);

void ellipticOccaRunTri2D(mesh2D *mesh);

void ellipticSetupTri2D(mesh2D *mesh, ogs_t **ogs, precon_t **precon, dfloat lambda);

void ellipticVolumeTri2D(mesh2D *mesh);

void ellipticSurfaceTri2D(mesh2D *mesh, dfloat time);

void ellipticUpdateTri2D(mesh2D *mesh, dfloat rka, dfloat rkb);

void ellipticErrorTri2D(mesh2D *mesh, dfloat time);

void ellipticParallelGatherScatter2D(mesh2D *mesh, ogs_t *ogs, occa::memory &o_v, occa::memory &o_gsv,
				     const char *type, const char *op);

precon_t *ellipticPreconditionerSetupTri2D(mesh2D *mesh, ogs_t *ogs, dfloat lambda);

void diagnostic(int N, occa::memory &o_x, const char *message);

void ellipticCoarsePreconditionerTri2D(mesh_t *mesh, precon_t *precon, dfloat *x, dfloat *b);

void ellipticCoarsePreconditionerSetupTri2D(mesh_t *mesh, precon_t *precon, dfloat lambda);
