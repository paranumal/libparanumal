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
  
  ogs_t *ogsP, *ogsDg;

  occa::memory o_diagA;
  
} precon_t;

void ellipticRunQuad2D(mesh2D *mesh);

void ellipticOccaRunQuad2D(mesh2D *mesh);

void ellipticSetupQuad2D(mesh2D *mesh, ogs_t **ogs, precon_t **precon, dfloat lambda);

void ellipticVolumeQuad2D(mesh2D *mesh);

void ellipticSurfaceQuad2D(mesh2D *mesh, dfloat time);

void ellipticUpdateQuad2D(mesh2D *mesh, dfloat rka, dfloat rkb);

void ellipticErrorQuad2D(mesh2D *mesh, dfloat time);

void ellipticParallelGatherScatter2D(mesh2D *mesh, ogs_t *ogs, occa::memory &o_v, occa::memory &o_gsv,
				     const char *type, const char *op);

precon_t *ellipticPreconditionerSetupQuad2D(mesh2D *mesh, ogs_t *ogs, dfloat lambda);

void diagnostic(int N, occa::memory &o_x, const char *message);
