#ifndef AGMG_H
#define AGMG_H 1

#ifdef OCCA_VERSION_1_0
#include "occa/modes/opencl/utils.hpp"
#endif

#include "mesh.h"

#include "parAlmond.h"
#include "agmgLevel.h"
#include "agmgMatrices.h"
#include "vectorPrimitives.h"
#include <mpi.h>

#define AGMGBDIM 32 //block size
#define SIMDWIDTH 32 //width of simd blocks
#define MAX_LEVELS 100
#define GPU_CPU_SWITCH_SIZE 0 //host-device switch threshold

#define RDIMX 32
#define RDIMY 8
#define RLOAD 1


void agmgSetup(parAlmond_t *parAlmond, csr *A, dfloat *nullA, hlong *globalRowStarts, setupAide options);
void parAlmondReport(parAlmond_t *parAlmond);
void buildAlmondKernels(parAlmond_t *parAlmond);

void kcycle(parAlmond_t *parAlmond, int k);
void device_kcycle(parAlmond_t *parAlmond, int k);

void vcycle(parAlmond_t *parAlmond, int k);
void device_vcycle(parAlmond_t *parAlmond, int k);

void pgmres(parAlmond_t *parAlmond, int maxIt, dfloat tol);
void device_pgmres(parAlmond_t *parAlmond, int maxIt, dfloat tol);

void pcg(parAlmond_t *parAlmond, int maxIt, dfloat tol);
void device_pcg(parAlmond_t *parAlmond, int maxIt, dfloat tol);

namespace agmg {
  extern int rank;
  extern int size;
  extern MPI_Comm comm;
};
#endif
