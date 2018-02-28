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
  occa::memory o_invDegreeP;

  occa::memory o_oasForwardDg;
  occa::memory o_oasBackDg;
  occa::memory o_oasDiagInvOpDg;
  occa::memory o_invDegreeDGP;
  
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

void ellipticSetupHex3D(mesh3D *mesh, occa::kernelInfo &kernelInfo);

void ellipticVolumeHex3D(mesh3D *mesh);

void ellipticSurfaceHex3D(mesh3D *mesh, dfloat time);

void ellipticUpdateHex3D(mesh3D *mesh, dfloat rka, dfloat rkb);

void ellipticErrorHex3D(mesh3D *mesh, dfloat time);

void ellipticParallelGatherScatterHex3D(mesh3D *mesh, ogs_t *ogs, occa::memory &o_v, occa::memory &o_gsv,
				     const char *type, const char *op);

precon_t *ellipticPreconditionerSetupHex3D(mesh3D *mesh, ogs_t *ogs, dfloat lambda, const char *options);

void ellipticCoarsePreconditionerHex3D(mesh_t *mesh, precon_t *precon, dfloat *x, dfloat *b);

void ellipticCoarsePreconditionerSetupHex3D(mesh_t *mesh, precon_t *precon, ogs_t *ogs, dfloat lambda, const char *options);

typedef struct {

  mesh_t *mesh;

  precon_t *precon;

  ogs_t *ogs;

  ogs_t *ogsDg;

  iint Nblock;
  
  occa::memory o_p; // search direction
  occa::memory o_z; // preconditioner solution
  occa::memory o_zP; // extended OAS preconditioner patch solution
  occa::memory o_Ax; // A*initial guess
  occa::memory o_Ap; // A*search direction
  occa::memory o_tmp; // temporary
  occa::memory o_grad; // temporary gradient storage (part of A*)
  occa::memory o_rtmp;
  occa::memory o_invDegree;
  
  dfloat *sendBuffer, *recvBuffer;

  // HOST shadow copies
  dfloat *Ax, *p, *r, *z, *zP, *Ap, *tmp, *grad;
  
}solver_t;

// block size for reduction (hard coded)
#define blockSize 256 

void ellipticMatrixFreeAx(void **args, occa::memory o_q, occa::memory o_Aq, const char* options);

int ellipticSolveHex3D(solver_t *solver, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const char *options);

solver_t *ellipticSolveSetupHex3D(mesh_t *mesh, dfloat lambda, occa::kernelInfo &kernelInfo, const char *options);