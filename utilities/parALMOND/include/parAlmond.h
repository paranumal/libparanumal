#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <occa.hpp>

#define AGMGBDIM 32 //block size
#define SIMDWIDTH 32 //width of simd blocks
#define MAX_LEVELS 100

#if 0
#define iint int
#define dfloat float
#define MPI_IINT MPI_INT
#define MPI_DFLOAT MPI_FLOAT
#define iintFormat "%d"
#define dfloatFormat "%f"
#define dfloatString "float"
#define iintString "int"
#else
#define iint int
#define dfloat double
#define MPI_IINT MPI_INT
#define MPI_DFLOAT MPI_DOUBLE
#define iintFormat "%d"
#define dfloatFormat "%lf"
#define dfloatString "double"
#define iintString "int"
#endif

typedef enum {PCG=0,GMRES=1}KrylovType;
typedef enum {JACOBI=0,DAMPED_JACOBI=1,CHEBYSHEV=2}SmoothType;

#include "parAlmondStructs.h"

typedef struct parAlmond {
	agmgLevel *levels;
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
  occa::kernel copyKernel; 
  occa::kernel scaleVectorKernel;
  occa::kernel partialInnerProdKernel;
  occa::kernel vectorAddKernel;
  occa::kernel vectorAddKernel2; 
  occa::kernel dotStarKernel; 
  occa::kernel simpleDotStarKernel;
  occa::kernel haloExtract; 

} almond_t;

#include "agmgLevel.h"
#include "coo.h"
#include "ell.h"
#include "hyb.h"
#include "dcsr.h"
#include "csr.h"
#include "vectorPrimitives.h"





void setup(csr *A, dfloat *nullA, iint *globalRowStarts);
void sync_setup_on_device(almond_t *almond, occa::device dev);
void buildAlmondKernels(almond_t *almond);

void solve(dfloat *rhs, dfloat* x, almond_t *almond);
void solve(occa::memory o_rhs, occa::memory o_x, almond_t *almond);

void kcycle(almond_t *almond, int k);
void device_kcycle(almond_t *almond, int k);

void vcycle(almond_t *almond, int k);
void device_vcycle(almond_t *almond, int k);

void gmres(almond_t *almond, csr *A, dfloat *b, dfloat *x, iint maxIt, dfloat tol);
void gmres(almond_t *almond, hyb *A, occa::memory o_b, occa::memory o_x, iint maxIt, dfloat tol);

void pcg(almond_t *almond, csr *A, dfloat *b, dfloat *x, iint maxIt, dfloat tol);
void pcg(almond_t *almond, hyb *A, occa::memory o_b, occa::memory o_x, iint maxIt, dfloat tol);