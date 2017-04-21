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
  void *xxt, *amg, *almond;

  occa::memory o_coarseInvDegree;
  occa::memory o_ztmp;

  iint coarseNp;
  iint coarseTotal;
  iint *coarseOffsets;
  dfloat *B, *tmp2;
  occa::memory *o_B, o_tmp2;
  void *xxt2;
  void *parAlmond;
  
} precon_t;

void ellipticRunHex3D(mesh3D *mesh);

void ellipticOccaRunHex3D(mesh3D *mesh);

void ellipticSetupHex3D(mesh3D *mesh, ogs_t **ogs, precon_t **precon, dfloat lambda, const char *options);

void ellipticVolumeHex3D(mesh3D *mesh);

void ellipticSurfaceHex3D(mesh3D *mesh, dfloat time);

void ellipticUpdateHex3D(mesh3D *mesh, dfloat rka, dfloat rkb);

void ellipticErrorHex3D(mesh3D *mesh, dfloat time);

void ellipticParallelGatherScatter3D(mesh3D *mesh, ogs_t *ogs, occa::memory &o_v, occa::memory &o_gsv,
				     const char *type, const char *op);

precon_t *ellipticPreconditionerSetupHex3D(mesh3D *mesh, ogs_t *ogs, dfloat lambda);

void ellipticCoarsePreconditionerHex3D(mesh_t *mesh, precon_t *precon, dfloat *x, dfloat *b);

void ellipticCoarsePreconditionerSetupHex3D(mesh_t *mesh, precon_t *precon, ogs_t *ogs, dfloat lambda, const char *options);

void ellipticOperator3D(mesh3D *mesh, dfloat *sendBuffer, dfloat *recvBuffer,
      ogs_t *ogs, dfloat lambda,
      occa::memory &o_q, occa::memory &o_gradq, occa::memory &o_Aq, const char *options);