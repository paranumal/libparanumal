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

  char *type;

  iint Nblock;

  dfloat tau;

  // HOST shadow copies
  dfloat *Ax, *p, *r, *z, *Ap, *tmp, *grad;

  dfloat *sendBuffer, *recvBuffer;

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


  occa::kernel AxKernel;
  occa::kernel rhsBCKernel;
  occa::kernel innerProductKernel;
  occa::kernel weightedInnerProduct1Kernel;
  occa::kernel weightedInnerProduct2Kernel;
  occa::kernel scaledAddKernel;
  occa::kernel dotMultiplyKernel;
  occa::kernel dotDivideKernel;

  occa::kernel gradientKernel;
  occa::kernel ipdgKernel;
  occa::kernel rhsBCIpdgKernel;

}solver_t;

void ellipticOccaRunTri2D(mesh2D *mesh);

void ellipticSetupTri2D(mesh2D *mesh, occa::kernelInfo &kernelInfo);

void ellipticParallelGatherScatterTri2D(mesh2D *mesh, ogs_t *ogs, occa::memory &o_v, occa::memory &o_gsv,
					const char *type, const char *op);

void ellipticPreconditionerSetupTri2D(solver_t *solver, ogs_t *ogs, dfloat tau, dfloat lambda, iint *BCType, const char *options, const char *parAlmondOptions);

void diagnostic(int N, occa::memory &o_x, const char *message);

void ellipticMatrixFreeAx(void **args, occa::memory o_q, occa::memory o_Aq, const char* options);

void ellipticPreconditioner2D(solver_t *solver, dfloat lambda, occa::memory &o_r,occa::memory &o_z,const char *options);

int ellipticSolveTri2D(solver_t *solver, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const char *options);

solver_t *ellipticSolveSetupTri2D(mesh_t *mesh, dfloat tau, dfloat lambda, iint *BCType, occa::kernelInfo &kernelInfo, const char *options, const char *parAlmondOptions);

solver_t *ellipticBuildMultigridLevelTri2D(solver_t *baseSolver, int* levelDegrees, int n, const char *options);

//smoother setups
void ellipticSetupSmootherOverlappingPatchIpdg(solver_t *solver, precon_t *precon, agmgLevel *level, dfloat tau, dfloat lambda, int *BCType, const char *options);
void ellipticSetupSmootherExactFullPatchIpdg  (solver_t *solver, precon_t *precon, agmgLevel *level, dfloat tau, dfloat lambda, int *BCType, const char *options);
void ellipticSetupSmootherApproxFullPatchIpdg (solver_t *solver, precon_t *precon, agmgLevel *level, dfloat tau, dfloat lambda, int *BCType, const char *options);
void ellipticSetupSmootherDampedJacobiIpdg    (solver_t *solver, precon_t *precon, agmgLevel *level, dfloat tau, dfloat lambda, int* BCType, const char *options);

void ellipticMultiGridSetupTri2D(solver_t *solver, precon_t* precon, dfloat tau, dfloat lambda, iint *BCType, const char *options, const char *parAlmondOptions);
void ellipticSetupSmootherTri2D(solver_t *solver, precon_t *precon,
                                dfloat tau, dfloat lambda, int* BCType,
                                const char *options);
dfloat maxEigSmoothAx(solver_t* solver, agmgLevel *level);