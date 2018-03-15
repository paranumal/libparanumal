#ifndef ELLIPTICHEX3D_H
#define ELLIPTICHEX3D_H 1

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh3D.h"
#include "parAlmond.h"
#include "ellipticPreconHex3D.h"

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

  int Nblock;

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

  //C0-FEM mask data
  int *mapB;      // boundary flag of face nodes
  dlong Nmasked;
  dlong *maskIds;

  occa::memory o_maskIds;
  occa::memory o_mapB;

  dfloat *sendBuffer, *recvBuffer;
  dfloat *gradSendBuffer, *gradRecvBuffer;

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
  dlong NglobalGatherElements;
  occa::memory o_globalGatherElementList;

  // list of elements that are not needed for global gather-scatter
  dlong NlocalGatherElements;
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

}solver_t;


void ellipticSetupHex3D(mesh3D *mesh, occa::kernelInfo &kernelInfo);

void ellipticParallelGatherScatterHex3D(mesh3D *mesh, ogs_t *ogs, occa::memory &o_v, const char *type, const char *op);

void ellipticPreconditionerSetupHex3D(solver_t *solver, ogs_t *ogs, dfloat tau, dfloat lambda, int *BCType, const char *options, const char *parAlmondOptions);

void diagnostic(int N, occa::memory &o_x, const char *message);

void ellipticPreconditioner3D(solver_t *solver, dfloat lambda, occa::memory &o_r,occa::memory &o_z,const char *options);

int ellipticSolveHex3D(solver_t *solver, dfloat lambda, dfloat tol, occa::memory &o_r, occa::memory &o_x, const char *options);

solver_t *ellipticSolveSetupHex3D(mesh_t *mesh, dfloat tau, dfloat lambda, int *BCType, occa::kernelInfo &kernelInfo, const char *options, const char *parAlmondOptions);

solver_t *ellipticBuildMultigridLevelHex3D(solver_t *baseSolver, int Nc, int Nf, int *BCType, const char *options);

void ellipticStartHaloExchange3D(solver_t *solver, occa::memory &o_q, int Nentries, dfloat *sendBuffer, dfloat *recvBuffer);

void ellipticInterimHaloExchange3D(solver_t *solver, occa::memory &o_q, int Nentries, dfloat *sendBuffer, dfloat *recvBuffer);

void ellipticEndHaloExchange3D(solver_t *solver, occa::memory &o_q, int Nentries, dfloat *recvBuffer);

void ellipticParallelGatherScatterSetup(solver_t *solver,const char *options);

//Linear solvers
int pcg      (solver_t* solver, const char* options, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const dfloat tol, const int MAXIT);

void ellipticScaledAdd(solver_t *solver, dfloat alpha, occa::memory &o_a, dfloat beta, occa::memory &o_b);
dfloat ellipticWeightedInnerProduct(solver_t *solver,
            occa::memory &o_w,
            occa::memory &o_a,
            occa::memory &o_b,
            const char *options);
void ellipticOperator3D(solver_t *solver, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq, const char *options);


void ellipticBuildIpdgHex3D(mesh3D *mesh, dfloat tau, dfloat lambda, int *BCType, nonZero_t **A,
                              dlong *nnzA, hlong *globalStarts, const char *options);

void ellipticBuildContinuousHex3D(solver_t *solver, dfloat lambda, nonZero_t **A, dlong *nnz,
                              ogs_t **ogs, hlong *globalStarts, const char* options);

void ellipticBuildJacobiHex3D(solver_t *solver, dfloat tau, dfloat lambda,
                                   int *BCType, dfloat **invDiagA,
                                   const char *options);

void ellipticBuildLocalPatchesHex3D(solver_t *solver, dfloat tau, dfloat lambda, int *BCType, dfloat rateTolerance,
                                   dlong *Npataches, dlong **patchesIndex, dfloat **patchesInvA,
                                   const char *options);

//smoother setups
void ellipticSetupSmootherOverlappingPatch(solver_t *solver, precon_t *precon, agmgLevel *level, dfloat tau, dfloat lambda, int *BCType, const char *options);
void ellipticSetupSmootherDampedJacobi    (solver_t *solver, precon_t *precon, agmgLevel *level, dfloat tau, dfloat lambda, int* BCType, const char *options);
void ellipticSetupSmootherLocalPatch(solver_t *solver, precon_t *precon, agmgLevel *level, dfloat tau, dfloat lambda, int *BCType, dfloat rateTolerance, const char *options);


void ellipticMultiGridSetupHex3D(solver_t *solver, precon_t* precon, dfloat tau, dfloat lambda, int *BCType, const char *options, const char *parAlmondOptions);
void ellipticSetupSmootherHex3D(solver_t *solver, precon_t *precon,
                                dfloat tau, dfloat lambda, int* BCType,
                                const char *options);

void matrixInverse(int N, dfloat *A);
dfloat maxEigSmoothAx(solver_t* solver, agmgLevel *level, const char* options);

void ellipticSEMFEMSetupHex3D(solver_t *solver, precon_t* precon,
                              dfloat tau, dfloat lambda, int *BCType,
                              const char *options, const char *parAlmondOptions);

#endif