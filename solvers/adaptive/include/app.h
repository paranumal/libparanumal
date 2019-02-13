#ifndef APP_H
#define APP_H 1

#include "adaptive.h"

typedef struct app
{
  occa::device device;

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
} app_t;

app_t *app_new(setupAide &options, MPI_Comm comm);
void app_free(app_t *app);

#endif
