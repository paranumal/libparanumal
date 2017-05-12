#ifndef PARALMOND_H 
#define PARALMOND_H 1

#include "mesh.h"
#include <mpi.h>

#define AGMGBDIM 32 //block size
#define SIMDWIDTH 32 //width of simd blocks
#define MAX_LEVELS 100

typedef enum {PCG=0,GMRES=1}KrylovType;
typedef enum {JACOBI=0,DAMPED_JACOBI=1,CHEBYSHEV=2}SmoothType;

#include "structs.h"

typedef struct {
	agmgLevel **levels;
	int numLevels;

	KrylovType ktype;

	//Coarse solver
	void (*coarseSolve)(void*,void*,void*);
	void *ACoarse;
	iint coarseTotal;
	iint coarseOffset;

	occa::device device;

  occa::kernel ellAXPYKernel; 
  occa::kernel ellZeqAXPYKernel;
  occa::kernel ellJacobi1Kernel;
  occa::kernel cooAXKernel1; 
  occa::kernel cooAXKernel2; 
  occa::kernel dcsrAXPYKernel; 
  occa::kernel dcsrZeqAXPYKernel; 
  occa::kernel dcsrJacobiKernel; 
  occa::kernel copyKernel; 
  occa::kernel scaleVectorKernel;
  occa::kernel partialInnerProdKernel;
  occa::kernel vectorAddKernel;
  occa::kernel vectorAddKernel2; 
  occa::kernel dotStarKernel; 
  occa::kernel simpleDotStarKernel;
  occa::kernel haloExtract;
  occa::kernel agg_interpolateKernel;

} almond_t;

#include "agmgLevel.h"
#include "agmgMatrices.h"
#include "vectorPrimitives.h"


almond_t *setup(csr *A, dfloat *nullA, iint *globalRowStarts);
void sync_setup_on_device(almond_t *almond, occa::device dev);
void buildAlmondKernels(almond_t *almond);

void solve(almond_t *almond, dfloat *rhs, dfloat* x);
void solve(almond_t *almond, occa::memory o_rhs, occa::memory o_x);

void kcycle(almond_t *almond, int k);
void device_kcycle(almond_t *almond, int k);

void vcycle(almond_t *almond, int k);
void device_vcycle(almond_t *almond, int k);

void gmres(almond_t *almond, csr *A, dfloat *b, dfloat *x, iint maxIt, dfloat tol);
void gmres(almond_t *almond, hyb *A, occa::memory o_b, occa::memory o_x, iint maxIt, dfloat tol);

void pcg(almond_t *almond, csr *A, dfloat *b, dfloat *x, iint maxIt, dfloat tol);
void pcg(almond_t *almond, hyb *A, occa::memory o_b, occa::memory o_x, iint maxIt, dfloat tol);

#endif