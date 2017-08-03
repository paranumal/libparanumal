#ifndef AGMG_H
#define AGMG_H 1

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


void agmgSetup(parAlmond_t *parAlmond, int lev, csr *A, dfloat *nullA, iint *globalRowStarts, const char* options);
void parAlmondReport(parAlmond_t *parAlmond);
void buildAlmondKernels(parAlmond_t *parAlmond);

void kcycle(parAlmond_t *parAlmond, int k);
void device_kcycle(parAlmond_t *parAlmond, int k);

void vcycle(parAlmond_t *parAlmond, int k);
void device_vcycle(parAlmond_t *parAlmond, int k);

void gmres(parAlmond_t *parAlmond, csr *A, dfloat *b, dfloat *x, iint maxIt, dfloat tol);
void gmres(parAlmond_t *parAlmond, hyb *A, occa::memory o_b, occa::memory o_x, iint maxIt, dfloat tol);

void pcg(parAlmond_t *parAlmond, iint maxIt, dfloat tol);
void device_pcg(parAlmond_t *parAlmond, iint maxIt, dfloat tol);

#endif