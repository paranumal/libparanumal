#ifndef ELLIPTICTRI2D_H
#define ELLIPTICTRI2D_H 1

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh2D.h"
#include "parAlmond.h"
#include "ellipticPreconTri2D.h"

// block size for reduction (hard coded)
#define blockSize 256

typedef struct {

  mesh_t *mesh;

  precon_t *precon;

  ogs_t *ogs;

  ogs_t *ogsDg;

  // C0 halo gather-scatter info
  ogs_t *halo;

  // C0 nonhalo gather-scatter info
  ogs_t *nonHalo;

  char *type;

  iint Nblock;

  dfloat tau;

  int *BCType;

  bool allNeumann;
  dfloat allNeumannPenalty;
  dfloat allNeumannScale;

  // HOST shadow copies
  dfloat *Ax, *p, *r, *z, *Ap, *tmp, *grad;
  dfloat *invDegree;

  dfloat *Ry, *R; //multigrid restriction matrix

  int *EToB;
  dfloat *sendBuffer, *recvBuffer;
  dfloat *gradSendBuffer, *gradRecvBuffer;

  //GMRES variables
  int GMRESrestartFreq;
  dfloat *HH;
  occa::memory *o_V;  

  occa::stream defaultStream;
  occa::stream dataStream;

  occa::memory o_p; // search direction
  occa::memory o_z; // preconditioner solution
  occa::memory o_res;
  occa::memory o_Sres;
  occa::memory o_Ax; // A*initial guess
  occa::memory o_Ap; // A*search direction
  occa::memory o_tmp; // temporary
  occa::memory o_grad; // temporary gradient storage (part of A*)
  occa::memory o_rtmp;
  occa::memory o_invDegree;
  occa::memory o_EToB;
  occa::memory o_R;
  occa::memory o_Ry;

  // list of elements that are needed for global gather-scatter
  iint NglobalGatherElements;
  occa::memory o_globalGatherElementList;

  // list of elements that are not needed for global gather-scatter
  iint NlocalGatherElements;
  occa::memory o_localGatherElementList;

  occa::kernel AxKernel;
  occa::kernel partialAxKernel;
  occa::kernel rhsBCKernel;
  occa::kernel addBCKernel;
  occa::kernel innerProductKernel;
  occa::kernel weightedInnerProduct1Kernel;
  occa::kernel weightedInnerProduct2Kernel;
  occa::kernel scaledAddKernel;
  occa::kernel dotMultiplyKernel;
  occa::kernel dotDivideKernel;

  occa::kernel gradientKernel;
  occa::kernel ipdgKernel;
  occa::kernel partialGradientKernel;
  occa::kernel partialIpdgKernel;
  occa::kernel rhsBCIpdgKernel;

  occa::kernel BRGradientVolumeKernel;
  occa::kernel BRGradientSurfaceKernel;
  occa::kernel BRDivergenceVolumeKernel;
  occa::kernel BRDivergenceSurfaceKernel;

}solver_t;

void ellipticSetupTri2D(mesh2D *mesh, occa::kernelInfo &kernelInfo, const char *options);

void ellipticParallelGatherScatterTri2D(mesh2D *mesh, ogs_t *ogs, occa::memory &o_v, occa::memory &o_gsv,
					const char *type, const char *op);

void ellipticPreconditionerSetupTri2D(solver_t *solver, ogs_t *ogs, dfloat tau, dfloat lambda, iint *BCType, const char *options, const char *parAlmondOptions);

void diagnostic(int N, occa::memory &o_x, const char *message);

void ellipticPreconditioner2D(solver_t *solver, dfloat lambda, occa::memory &o_r,occa::memory &o_z,const char *options);

int ellipticSolveTri2D(solver_t *solver, dfloat lambda, dfloat tol, occa::memory &o_r, occa::memory &o_x, const char *options);

solver_t *ellipticSolveSetupTri2D(mesh_t *mesh, dfloat tau, dfloat lambda, iint *BCType, occa::kernelInfo &kernelInfo, const char *options, const char *parAlmondOptions);

solver_t *ellipticBuildMultigridLevelTri2D(solver_t *baseSolver, int Nc, int Nf, int* BCType, const char *options);

void ellipticStartHaloExchange2D(solver_t *solver, occa::memory &o_q, int Nentries, dfloat *sendBuffer, dfloat *recvBuffer);

void ellipticInterimHaloExchange2D(solver_t *solver, occa::memory &o_q, int Nentries, dfloat *sendBuffer, dfloat *recvBuffer);

void ellipticEndHaloExchange2D(solver_t *solver, occa::memory &o_q, int Nentries, dfloat *recvBuffer);

void ellipticParallelGatherScatterSetup(solver_t *solver,const char *options);

//Linear solvers
int pcg      (solver_t* solver, const char* options, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const dfloat tol, const int MAXIT);
int pbicgstab(solver_t* solver, const char* options, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const dfloat tol, const int MAXIT);
int pgmresm  (solver_t* solver, const char* options, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const dfloat tol, const int MAXIT);


dfloat ellipticScaledAdd(solver_t *solver, dfloat alpha, occa::memory &o_a, dfloat beta, occa::memory &o_b);
dfloat ellipticWeightedInnerProduct(solver_t *solver,
            occa::memory &o_w,
            occa::memory &o_a,
            occa::memory &o_b,
            const char *options);
void ellipticOperator2D(solver_t *solver, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq, const char *options);


void ellipticBuildJacobiTri2D(solver_t *solver, mesh2D *mesh, iint basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda,
                                   iint *BCType, dfloat **invDiagA,
                                   const char *options);

void ellipticBuildFullPatchesTri2D(solver_t *solver, mesh2D *mesh, iint basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda, iint *BCType, dfloat rateTolerance,
                                   iint *Npataches, iint **patchesIndex, dfloat **patchesInvA,
                                   const char *options);

void ellipticBuildFacePatchesTri2D(solver_t *solver, mesh2D *mesh, iint basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda, iint *BCType, dfloat rateTolerance,
                                   iint *Npataches, iint **patchesIndex, dfloat **patchesInvA,
                                   const char *options);

void ellipticBuildLocalPatchesTri2D(solver_t *solver, mesh2D *mesh, iint basisNp, dfloat *basis,
                                   dfloat tau, dfloat lambda, iint *BCType, dfloat rateTolerance,
                                   iint *Npataches, iint **patchesIndex, dfloat **patchesInvA,
                                   const char *options);

//smoother setups
void ellipticSetupSmootherOverlappingPatch(solver_t *solver, precon_t *precon, agmgLevel *level, dfloat tau, dfloat lambda, int *BCType, const char *options);
void ellipticSetupSmootherDampedJacobi    (solver_t *solver, precon_t *precon, agmgLevel *level, dfloat tau, dfloat lambda, int* BCType, const char *options);
void ellipticSetupSmootherFullPatch (solver_t *solver, precon_t *precon, agmgLevel *level, dfloat tau, dfloat lambda, int *BCType, dfloat rateTolerance, const char *options);
void ellipticSetupSmootherFacePatch (solver_t *solver, precon_t *precon, agmgLevel *level, dfloat tau, dfloat lambda, int *BCType, dfloat rateTolerance, const char *options);
void ellipticSetupSmootherLocalPatch(solver_t *solver, precon_t *precon, agmgLevel *level, dfloat tau, dfloat lambda, int *BCType, dfloat rateTolerance, const char *options);


void ellipticMultiGridSetupTri2D(solver_t *solver, precon_t* precon, dfloat tau, dfloat lambda, iint *BCType, const char *options, const char *parAlmondOptions);
void ellipticSetupSmootherTri2D(solver_t *solver, precon_t *precon,
                                dfloat tau, dfloat lambda, int* BCType,
                                const char *options);

void matrixInverse(int N, dfloat *A);
dfloat maxEigSmoothAx(solver_t* solver, agmgLevel *level, const char* options);

void ellipticSEMFEMSetupTri2D(solver_t *solver, precon_t* precon,
                              dfloat tau, dfloat lambda, iint *BCType,
                              const char *options, const char *parAlmondOptions);

#endif

//KS: load sparse stiffness ops
void loadElementStiffnessMatricesTri2D(mesh2D *mesh, const char *options, int N);
void buildElementStiffnessMatricesTri2D(mesh2D *mesh, const char *options, int N);
