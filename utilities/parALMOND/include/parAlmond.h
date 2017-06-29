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

  mesh_t *mesh;
  hgs_t *hgs;
  const char* options;

  //Matrix Free args 
  void (*MatFreeAx)(void **args, occa::memory o_q, occa::memory o_Aq,const char* options);
  void **MatFreeArgs;

	//Coarse xxt solver
	void *ExactSolve;
	iint coarseTotal;
	iint coarseOffset;
  dfloat *xCoarse, *rhsCoarse;

	occa::device device;

  occa::memory o_x;
  occa::memory o_Ax;
  occa::memory o_rho;

  occa::kernel ellAXPYKernel; 
  occa::kernel ellZeqAXPYKernel;
  occa::kernel ellJacobiKernel;
  occa::kernel cooAXKernel; 
  occa::kernel scaleVectorKernel;
  occa::kernel vectorAddKernel;
  occa::kernel vectorAddKernel2; 
  occa::kernel setVectorKernel;
  occa::kernel dotStarKernel; 
  occa::kernel simpleDotStarKernel;
  occa::kernel haloExtract;
  occa::kernel agg_interpolateKernel;
  occa::kernel innerProdKernel;
  occa::kernel vectorAddInnerProdKernel;
  occa::kernel kcycleCombinedOp1Kernel;
  occa::kernel kcycleCombinedOp2Kernel;

} parAlmond_t;

#include "agmgLevel.h"
#include "agmgMatrices.h"
#include "vectorPrimitives.h"


parAlmond_t *parAlmondInit(mesh_t *mesh, const char* options);

void *parAlmondAgmgSetup(parAlmond_t *parAlmond,
       iint* rowStarts,
       iint  nnz,
       iint* Ai,
       iint* Aj,
       dfloat* Avals,
       hgs_t *hgs,
       const char* options);

void parAlmondPrecon(occa::memory o_x,
    void* ALMOND,
    occa::memory o_rhs);

int parAlmondFree(void* A);

void parAlmondSetMatFreeAX(void* A, void (*MatFreeAx)(void **args, occa::memory o_q, occa::memory o_Aq,const char* options),
                        void **args);

void agmgSetup(parAlmond_t *parAlmond, csr *A, dfloat *nullA, iint *globalRowStarts, const char* options);
void sync_setup_on_device(parAlmond_t *parAlmond);
void parAlmondReport(parAlmond_t *parAlmond);
void buildAlmondKernels(parAlmond_t *parAlmond);

void parAlmondMatrixFreeAX(parAlmond_t *parAlmond, occa::memory &o_x, occa::memory &o_Ax);

void kcycle(parAlmond_t *parAlmond, int k);
void device_kcycle(parAlmond_t *parAlmond, int k);

void vcycle(parAlmond_t *parAlmond, int k);
void device_vcycle(parAlmond_t *parAlmond, int k);

void gmres(parAlmond_t *parAlmond, csr *A, dfloat *b, dfloat *x, iint maxIt, dfloat tol);
void gmres(parAlmond_t *parAlmond, hyb *A, occa::memory o_b, occa::memory o_x, iint maxIt, dfloat tol);

void pcg(parAlmond_t *parAlmond, iint maxIt, dfloat tol);
void device_pcg(parAlmond_t *parAlmond, iint maxIt, dfloat tol);

#endif