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

// block size for reduction (hard coded)
#define blockSize 256

typedef struct {

  int dim;
  int elementType; // number of edges (3=tri, 4=quad, 6=tet, 12=hex)

  mesh_t *mesh;

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
  occa::memory h_sendBuffer, h_recvBuffer;

  occa::memory o_gradSendBuffer, o_gradRecvBuffer;
  occa::memory h_gradSendBuffer, h_gradRecvBuffer;

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

  occa::memory o_ggeoNoJW; // 

  occa::kernel AxKernel;
  occa::kernel partialAxKernel;
  occa::kernel partialAxKernel2;
  occa::kernel partialFloatAxKernel;
  occa::kernel partialCubatureAxKernel;
  
  occa::kernel rhsBCKernel;
  occa::kernel addBCKernel;
  occa::kernel innerProductKernel;
  occa::kernel weightedInnerProduct1Kernel;
  occa::kernel weightedInnerProduct2Kernel;
  occa::kernel scaledAddKernel;
  occa::kernel dotMultiplyKernel;
  occa::kernel dotMultiplyAddKernel;
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
  dlong           NblocksUpdatePCG;
  dfloat         *tmpNormr;
  occa::memory  o_tmpNormr;
  occa::kernel  updatePCGKernel;

  // combined NBPCG update steps
  dfloat *tmppdots;
  dfloat *tmprdotz;
  dfloat *tmpzdotz;
  dfloat *tmprdotr;
  occa::kernel update1NBPCGKernel;
  occa::kernel update2NBPCGKernel;
  occa::memory  o_tmppdots;
  occa::memory  o_tmprdotz;
  occa::memory  o_tmpzdotz;
  occa::memory  o_tmprdotr;
  //  occa::memory  o_s;
  //  occa::memory  o_S;
  //  occa::memory  o_Z;

  occa::memory *o_pcgWork;

  // combined NBFPCG update steps
  dfloat *tmpudotr;
  dfloat *tmpudots;
  dfloat *tmpudotw;
  occa::kernel update0NBFPCGKernel;
  occa::kernel update1NBFPCGKernel;
  occa::memory  o_tmpudotr;
  occa::memory  o_tmpudots;
  occa::memory  o_tmpudotw;

  hlong NelementsGlobal;
  dfloat nullProjectWeightGlobal;
  
}elliptic_t;

#include "ellipticMultiGrid.h"

elliptic_t *ellipticSetup(mesh2D *mesh, dfloat lambda, occa::properties &kernelInfo, setupAide options);

void ellipticPreconditioner(elliptic_t *elliptic, dfloat lambda, occa::memory &o_r, occa::memory &o_z);
void ellipticPreconditionerSetup(elliptic_t *elliptic, ogs_t *ogs, dfloat lambda,
				 occa::properties &kernelInfo);

int  ellipticSolve(elliptic_t *elliptic, dfloat lambda, dfloat tol, occa::memory &o_r, occa::memory &o_x);

void ellipticSolveSetup(elliptic_t *elliptic, dfloat lambda, occa::properties &kernelInfo);


void ellipticStartHaloExchange(elliptic_t *elliptic, occa::memory &o_q, int Nentries, dfloat *sendBuffer, dfloat *recvBuffer);
void ellipticInterimHaloExchange(elliptic_t *elliptic, occa::memory &o_q, int Nentries, dfloat *sendBuffer, dfloat *recvBuffer);
void ellipticEndHaloExchange(elliptic_t *elliptic, occa::memory &o_q, int Nentries, dfloat *recvBuffer);

//Linear solvers
int pcg      (elliptic_t* elliptic, dfloat lambda, occa::memory &o_r, occa::memory &o_x, const dfloat tol, const int MAXIT);



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

void ellipticSerialPartialAxHexKernel3D(const int Nq,
					const hlong Nelements,
					const occa::memory &o_elementList,
					const occa::memory &o_ggeo,
					const occa::memory &o_Dmatrices,
					const occa::memory &o_Smatrices,
					const occa::memory &o_MM,
					const dfloat lambda,
					const occa::memory &o_q,
					occa::memory &o_Aq);

void ellipticSerialAxHexKernel3D(const int Nq,
				 const hlong Nelements,
				 const occa::memory &o_ggeo,
				 const occa::memory &o_Dmatrices,
				 const occa::memory &o_Smatrices,
				 const occa::memory &o_MM,
				 const dfloat lambda,
				 const occa::memory &o_q,
				 occa::memory &o_Aq,
				 const occa::memory &o_ggeoNoJW);

void ellipticBuildOneRing(elliptic_t *elliptic, dfloat lambda, occa::properties &kernelInfo);
void ellipticOneRingDiagnostics(elliptic_t *elliptic, elliptic_t *elliptic1, dfloat lambda);

occa::properties ellipticKernelInfo(mesh_t *mesh);

void ellipticOasSetup(elliptic_t *elliptic, dfloat lambda,
		      occa::properties &kernelInfo);
void ellipticOasSolve(elliptic_t *elliptic, dfloat lambda,
		      occa::memory &o_r, occa::memory &o_x);

void ellipticOneRingExchange(elliptic_t *elliptic,
			     elliptic_t *elliptic1,
			     size_t Nbytes,       // message size per element
			     occa::memory &o_q,
			     occa::memory &o_qOneRing);


void ellipticNonBlockingUpdate1NBPCG(elliptic_t *elliptic,
				     occa::memory &o_z, occa::memory &o_Z, const dfloat beta,
				     occa::memory &o_p, occa::memory &o_s,
				     dfloat *localpdots, dfloat *globalpdots, MPI_Request *request);

void ellipticNonBlockingUpdate2NBPCG(elliptic_t *elliptic,
				     occa::memory &o_s, occa::memory &o_S, const dfloat alpha,
				     occa::memory &o_r, occa::memory &o_z,
				     dfloat *localdots, dfloat *globaldots, MPI_Request *request);

int nbpcg(elliptic_t* elliptic, dfloat lambda, 
	  occa::memory &o_r, occa::memory &o_x, 
	  const dfloat tol, const int MAXIT);

void ellipticNonBlockingUpdate0NBFPCG(elliptic_t *elliptic,
				      occa::memory &o_u, occa::memory &o_r, occa::memory &o_w,
				      dfloat *localdots, dfloat *globaldots, MPI_Request *request);

void ellipticNonBlockingUpdate1NBFPCG(elliptic_t *elliptic,
				      occa::memory &o_p, occa::memory &o_s, occa::memory &o_q, occa::memory &o_z,
				      const dfloat alpha,
				      occa::memory &o_x, occa::memory &o_r, occa::memory &o_u, occa::memory &o_w,
				      dfloat *localpdots, dfloat *globalpdots, MPI_Request *request);

int nbfpcg(elliptic_t* elliptic, dfloat lambda, 
	   occa::memory &o_r, occa::memory &o_x, 
	   const dfloat tol, const int MAXIT);

void ellipticZeroMean(elliptic_t *elliptic, occa::memory &o_q);

#endif

