/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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
