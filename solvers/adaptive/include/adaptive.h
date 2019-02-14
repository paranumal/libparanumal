#ifndef ADAPTIVE_H
#define ADAPTIVE_H 1

#include "headers.h"

typedef struct adaptive
{

  MPI_Comm comm;
  int rank;
  int size;
  
  occa::device device;

  int allNeumann; // One if all Neumann and lambda==0, else zero
  
  p4est_connectivity_t *conn;
  p4est_t *pxest;
  p4est_ghost_t *ghost;

  int brick_n[DIM];
  int brick_p[DIM];
  p4est_topidx_t *brick_TToC; // tree id to cartesian coordinates needed for
                              // periodic bricks

  level_t *lvl;

  // generic kernels [ e.g. PCG ]
  occa::kernel sumKernel;
  occa::kernel addScalarKernel;
  occa::kernel scaledAddKernel;
  occa::kernel dotMultiplyKernel;
  occa::kernel dotMultiplyAddKernel;
  occa::kernel dotDivideKernel;

  occa::kernel innerProductKernel;
  occa::kernel weightedInnerProduct1Kernel;
  occa::kernel weightedInnerProduct2Kernel;

  occa::kernel weightedNorm2Kernel;
  occa::kernel norm2Kernel;
  
  occa::kernel updatePCGKernel;
  occa::kernel update1NBPCGKernel;
  occa::kernel update2NBPCGKernel;
  occa::kernel update0NBFPCGKernel;
  occa::kernel update1NBFPCGKernel;

  setupAide options;
  
} adaptive_t;

adaptive_t *adaptive_new(setupAide &options, MPI_Comm comm);
void adaptive_free(adaptive_t *adaptive);

void adaptivePreconditioner(adaptive_t *adaptive, dfloat lambda,
                            occa::memory &o_r, occa::memory &o_z);

dfloat adaptiveWeightedNorm2(adaptive_t *adaptive, level_t *level,
			     occa::memory &o_w, occa::memory &o_a);

void adaptiveScaledAdd(adaptive_t *adaptive, level_t *level,
		       dfloat alpha, occa::memory &o_a, dfloat beta, occa::memory &o_b);

void adaptiveZeroMean(adaptive_t *adaptive, level_t *level, occa::memory &o_q);

#endif
