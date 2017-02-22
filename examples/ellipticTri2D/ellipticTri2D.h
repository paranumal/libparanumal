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

  occa::memory o_oasForwardDgT;
  occa::memory o_oasBackDgT;
  
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

void ellipticSetupTri2D(mesh2D *mesh, occa::kernelInfo &kernelInfo);

void ellipticVolumeTri2D(mesh2D *mesh);

void ellipticSurfaceTri2D(mesh2D *mesh, dfloat time);

void ellipticUpdateTri2D(mesh2D *mesh, dfloat rka, dfloat rkb);

void ellipticErrorTri2D(mesh2D *mesh, dfloat time);

void ellipticParallelGatherScatterTri2D(mesh2D *mesh, ogs_t *ogs, occa::memory &o_v, occa::memory &o_gsv,
					const char *type, const char *op);

precon_t *ellipticPreconditionerSetupTri2D(mesh2D *mesh, ogs_t *ogs, dfloat lambda);

void diagnostic(int N, occa::memory &o_x, const char *message);

void ellipticCoarsePreconditionerTri2D(mesh_t *mesh, precon_t *precon, dfloat *x, dfloat *b);

void ellipticCoarsePreconditionerSetupTri2D(mesh_t *mesh, precon_t *precon, dfloat lambda);


typedef struct {

  mesh_t *mesh;

  precon_t *precon;

  ogs_t *ogs;

  ogs_t *ogsDg;

  char *type;

  iint Nblock;
  
  occa::memory o_p; // search direction
  occa::memory o_z; // preconditioner solution
  occa::memory o_zP; // extended OAS preconditioner patch solution
  occa::memory o_Ax; // A*initial guess
  occa::memory o_Ap; // A*search direction
  occa::memory o_tmp; // temporary
  occa::memory o_grad; // temporary gradient storage (part of A*)
  occa::memory o_rtmp;
  
  dfloat *sendBuffer, *recvBuffer;

  // HOST shadow copies
  dfloat *Ax, *p, *r, *z, *zP, *Ap, *tmp, *grad;
  
}solver_t;

// block size for reduction (hard coded)
#define blockSize 256 


int ellipticSolveTri2D(solver_t *solver, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const char *options);

solver_t *ellipticSolveSetupTri2D(mesh_t *mesh, dfloat lambda, occa::kernelInfo &kernelInfo);
