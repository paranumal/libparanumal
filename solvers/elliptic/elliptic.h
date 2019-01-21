/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#ifndef ELLIPTIC_H
#define ELLIPTIC_H 1

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh2D.h"
#include "mesh3D.h"
#include "parAlmond.hpp"
#include "ellipticPrecon.h"
#include "timer.h"

// block size for reduction (hard coded)
#define blockSize 256

typedef struct {

  int DEBUG_ENABLE_REDUCTIONS; 
  int DEBUG_ENABLE_MEMCOPY; 
  int DEBUG_ENABLE_OGS; 
  int DEBUG_ENABLE_MPIREDUCE; 

  int dim;
  int elementType; // number of edges (3=tri, 4=quad, 6=tet, 12=hex)

  mesh_t *mesh;
  timer *profiler; 

  precon_t *precon;

  ogs_t *ogs;

  setupAide options;

  char *type;

  dlong Nblock;
  dlong Nblock2; // second reduction

  dfloat tau;

  int *BCType;

  bool allNeumann;
  dfloat allNeumannPenalty;
  dfloat allNeumannScale;

  // HOST shadow copies
  dfloat *x, *Ax, *p, *r, *z, *Ap, *tmp, *grad;
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
  occa::memory o_sendBuffer, o_recvBuffer;
  occa::memory o_gradSendBuffer, o_gradRecvBuffer;


  occa::stream defaultStream;
  occa::stream dataStream;

  occa::memory o_x;
  occa::memory o_r;
  occa::memory o_p; // search direction
  occa::memory o_z; // preconditioner solution
  occa::memory o_res;
  occa::memory o_Sres;
  occa::memory o_Ax; // A*initial guess
  occa::memory o_Ap; // A*search direction
  occa::memory o_tmp; // temporary
  occa::memory o_tmp2; // temporary (second reduction)
  occa::memory o_grad; // temporary gradient storage (part of A*)
  occa::memory o_rtmp;
  occa::memory o_invDegree;
  occa::memory o_EToB;
  occa::memory o_R;
  occa::memory o_Ry;

  occa::memory o_EXYZ; // element vertices for reconstructing geofacs (trilinear hexes only)
  occa::memory o_gllzw; // GLL nodes and weights

  occa::kernel AxKernel;
  occa::kernel partialAxKernel;
  occa::kernel partialFloatAxKernel;
  occa::kernel partialCubatureAxKernel;
  
  occa::kernel rhsBCKernel;
  occa::kernel addBCKernel;
  occa::kernel innerProductKernel;
  occa::kernel weightedInnerProduct1Kernel;
  occa::kernel weightedInnerProduct2Kernel;
  occa::kernel scaledAddKernel;
  occa::kernel dotMultiplyKernel;
  occa::kernel dotDivideKernel;

  occa::kernel weightedNorm2Kernel;
  occa::kernel norm2Kernel;

  occa::kernel gradientKernel;
  occa::kernel ipdgKernel;
  occa::kernel partialGradientKernel;
  occa::kernel partialIpdgKernel;
  occa::kernel rhsBCIpdgKernel;

  // combined PCG update step
  int             NthreadsUpdatePCG;
  hlong           NblocksUpdatePCG;
  dfloat         *tmpNormr;
  occa::memory  o_tmpNormr;
  occa::kernel  updatePCGKernel;
  
}elliptic_t;

#include "ellipticMultiGrid.h"

elliptic_t *ellipticSetup(mesh2D *mesh, dfloat lambda, occa::properties &kernelInfo, setupAide options);

void ellipticPreconditioner(elliptic_t *elliptic, dfloat lambda, occa::memory &o_r, occa::memory &o_z);
void ellipticPreconditionerSetup(elliptic_t *elliptic, ogs_t *ogs, dfloat lambda);

int  ellipticSolve(elliptic_t *elliptic, dfloat lambda, dfloat tol, occa::memory &o_r, occa::memory &o_x);
int  ellipticSolveTest(elliptic_t *elliptic, dfloat lambda, dfloat tol, occa::memory &o_r, occa::memory &o_x, int maxIter);
void ellipticSolveSetup(elliptic_t *elliptic, dfloat lambda, occa::properties &kernelInfo);


void ellipticStartHaloExchange(elliptic_t *elliptic, occa::memory &o_q, int Nentries, dfloat *sendBuffer, dfloat *recvBuffer);
void ellipticInterimHaloExchange(elliptic_t *elliptic, occa::memory &o_q, int Nentries, dfloat *sendBuffer, dfloat *recvBuffer);
void ellipticEndHaloExchange(elliptic_t *elliptic, occa::memory &o_q, int Nentries, dfloat *recvBuffer);

//Linear solvers
int pcg      (elliptic_t* elliptic, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const dfloat tol, const int MAXIT);
void pcgBP5   (elliptic_t* elliptic, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const int MAXIT);

void ellipticScaledAdd(elliptic_t *elliptic, dfloat alpha, occa::memory &o_a, dfloat beta, occa::memory &o_b);
dfloat ellipticWeightedInnerProduct(elliptic_t *elliptic, occa::memory &o_w, occa::memory &o_a, occa::memory &o_b);

dfloat ellipticCascadingWeightedInnerProduct(elliptic_t *elliptic, occa::memory &o_w, occa::memory &o_a, occa::memory &o_b);

void ellipticOperator(elliptic_t *elliptic, dfloat lambda, occa::memory &o_q, occa::memory &o_Aq, const char *precision);

dfloat ellipticWeightedNorm2(elliptic_t *elliptic, occa::memory &o_w, occa::memory &o_a);
void ellipticBuildIpdg(elliptic_t* elliptic, int basisNp, dfloat *basis, dfloat lambda,
                        nonZero_t **A, dlong *nnzA, hlong *globalStarts);

void ellipticBuildContinuous(elliptic_t* elliptic, dfloat lambda, nonZero_t **A,
                                  dlong *nnz, ogs_t **ogs, hlong *globalStarts);

void ellipticBuildJacobi(elliptic_t *elliptic, dfloat lambda, dfloat **invDiagA);

void ellipticBuildLocalPatches(elliptic_t *elliptic, dfloat lambda, dfloat rateTolerance,
                               dlong *Npataches, dlong **patchesIndex, dfloat **patchesInvA);

// //smoother setups
// void ellipticSetupSmoother(elliptic_t *elliptic, precon_t *precon, dfloat lambda);
// void ellipticSetupSmootherDampedJacobi    (elliptic_t *elliptic, precon_t *precon, agmgLevel *level, dfloat lambda);
// void ellipticSetupSmootherLocalPatch(elliptic_t *elliptic, precon_t *precon, agmgLevel *level, dfloat lambda, dfloat rateTolerance);

void ellipticMultiGridSetup(elliptic_t *elliptic, precon_t* precon, dfloat lambda);
elliptic_t *ellipticBuildMultigridLevel(elliptic_t *baseElliptic, int Nc, int Nf);

void ellipticSEMFEMSetup(elliptic_t *elliptic, precon_t* precon, dfloat lambda);

dfloat ellipticUpdatePCG(elliptic_t *elliptic, occa::memory &o_p, occa::memory &o_Ap, dfloat alpha,
			 occa::memory &o_x, occa::memory &o_r);

// dfloat maxEigSmoothAx(elliptic_t* elliptic, agmgLevel *level);

#define maxNthreads 256

extern "C"
{
  void ellipticPlotVTUHex3D(mesh3D *mesh, char *fileNameBase, int fld);
}

#endif

