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

#ifndef PARALMOND_H
#define PARALMOND_H 1

typedef struct csr_t {

  dlong Nrows;
  dlong Ncols;

  dlong NlocalCols;

  //local
  dlong diagNNZ;
  dlong   *diagRowStarts;
  dlong   *diagCols;
  dfloat *diagCoefs;

  //non-local
  dlong offdNNZ;
  dlong   *offdRowStarts;
  dlong   *offdCols;
  dfloat *offdCoefs;

  dfloat *diagInv;

  dfloat *null;

  //storage for smoothing
  dfloat *scratch;

  hlong   *colMap;

  // MPI halo exchange info
  dlong  NHalo;
  int  NrecvTotal;  // number of elements to be sent in halo exchange
  int  NsendTotal;
  dlong  totalHaloPairs;
  dlong *haloElementList; // sorted list of elements to be sent in halo exchange
  int *NsendPairs;      // number of elements worth of data to send
  int *NrecvPairs;      // number of elements worth of data to recv
  int  NsendMessages;   // number of messages to send
  int  NrecvMessages;   // number of messages to recv
  dfloat *sendBuffer;

  void *haloSendRequests;
  void *haloRecvRequests;

} csr;


typedef struct ell_t {

  dlong Nrows;
  dlong Ncols;
  int nnzPerRow;
  dlong strideLength;
  dlong actualNNZ;

  occa::memory o_cols;
  occa::memory o_coefs;

} ell;

typedef struct coo_t {

  dlong Nrows;
  dlong Ncols;
  dlong nnz;

  // device memory
  occa::memory o_offsets;
  occa::memory o_cols;
  occa::memory o_coefs;

} coo;

typedef struct hyb_t {

  dlong Nrows;
  dlong Ncols;

  dlong NlocalCols;

  coo *C;
  ell *E;

  occa::memory o_diagInv;

  occa::memory o_null;

  // MPI halo exchange info
  dlong  NHalo;
  hlong *colMap;
  int  NrecvTotal;  // number of elements to be sent in halo exchange
  int  NsendTotal;
  dlong *haloElementList; // sorted list of elements to be sent in halo exchange
  occa::memory o_haloElementList;
  int *NsendPairs;      // number of elements worth of data to send
  int *NrecvPairs;      // number of elements worth of data to recv
  int  NsendMessages;   // number of messages to send
  int  NrecvMessages;   // number of messages to recv
  dfloat   *sendBuffer;
  dfloat   *recvBuffer;
  occa::memory o_haloBuffer;

  void *haloSendRequests;
  void *haloRecvRequests;

} hyb;


typedef struct dcsr_t {

  dlong Nrows;
  dlong Ncols;

  dlong NlocalCols;

  //local
  dlong diagNNZ;
  occa::memory o_diagRows;
  occa::memory o_diagCols;
  occa::memory o_diagCoefs;

  //non-local
  dlong offdNNZ;
  occa::memory o_offdRows;
  occa::memory o_offdCols;
  occa::memory o_offdCoefs;

  // MPI halo exchange info
  dlong  NHalo;
  int  NrecvTotal;  // number of elements to be sent in halo exchange
  int  NsendTotal;
  dlong  totalHaloPairs;
  dlong *haloElementList; // sorted list of elements to be sent in halo exchange
  int *NsendPairs;      // number of elements worth of data to send
  int *NrecvPairs;      // number of elements worth of data to recv
  int  NsendMessages;   // number of messages to send
  int  NrecvMessages;   // number of messages to recv
  dfloat   *sendBuffer;
  dfloat   *recvBuffer;

  occa::memory o_haloElementList;
  occa::memory o_haloBuffer;

  void *haloSendRequests;
  void *haloRecvRequests;

} dcoo;

typedef enum {PCG=0,GMRES=1}KrylovType;
typedef enum {JACOBI=0,DAMPED_JACOBI=1,CHEBYSHEV=2}SmoothType;

typedef struct agmgLevel_t {
  dlong Nrows;
  dlong Ncols;

  hlong *globalRowStarts; //global partitioning of fine level
  hlong *globalAggStarts; //global partitioning of coarse level

  bool gatherLevel;
  bool weightedInnerProds;

  void **AxArgs;
  void **smoothArgs;
  void **smootherArgs;
  void **coarsenArgs;
  void **prolongateArgs;
  void **gatherArgs;
  void **scatterArgs;

  //operator call-backs
  void (*device_Ax)        (void **args, occa::memory &o_x, occa::memory &o_Ax);
  void (*device_smooth)    (void **args, occa::memory &o_r, occa::memory &o_x, bool x_is_zero);
  void (*device_smoother)  (void **args, occa::memory &o_r, occa::memory &o_Sr);
  void (*device_coarsen)   (void **args, occa::memory &o_x, occa::memory &o_Rx);
  void (*device_prolongate)(void **args, occa::memory &o_x, occa::memory &o_Px);
  void (*device_gather)    (void **args, occa::memory &o_x, occa::memory &o_Gx);
  void (*device_scatter)   (void **args, occa::memory &o_x, occa::memory &o_Sx);

  //host versions
  void (*Ax)        (void **args, dfloat *x, dfloat *Ax);
  void (*smooth)    (void **args, dfloat *r, dfloat *x, bool x_is_zero);
  void (*smoother)  (void **args, dfloat *r, dfloat *Sr);
  void (*coarsen)   (void **args, dfloat *x, dfloat *Rx);
  void (*prolongate)(void **args, dfloat *x, dfloat *Px);
  void (*gather)    (void **args, dfloat *x, dfloat *Gx);
  void (*scatter)   (void **args, dfloat *x, dfloat *Sx);

  //agmg operators
  csr *A;
  csr *P;
  csr *R;

  hyb  *deviceA;
  dcoo  *dcsrP;
  hyb  *deviceR;

  dfloat *rhs, *res, *x;

  dfloat *Srhs, *Sx; //scatter copies

  dfloat *ckp1, *vkp1, *wkp1;

  dfloat *weight;

  occa::memory o_rhs, o_res, o_x;
  occa::memory o_Srhs, o_Sx;
  occa::memory o_ckp1, o_vkp1, o_wkp1;

  occa::memory o_weight;

  dfloat *smoother_params;
  dfloat *smootherResidual;
  dfloat *smootherResidual2;
  dfloat *smootherUpdate;
  occa::memory o_smootherResidual;
  occa::memory o_smootherResidual2;
  occa::memory o_smootherUpdate;
  int ChebyshevIterations;

  dfloat threshold;
  dlong numAggregates;
  SmoothType stype;

} agmgLevel;

typedef struct {
  agmgLevel **levels;
  int numLevels;

  KrylovType ktype;

  setupAide options;

  //Matrix Free args
  void (*MatFreeAx)(void **args, occa::memory o_q, occa::memory o_Aq,const char* options);
  void **MatFreeArgs;

  //Coarse solver
  void *ExactSolve;
  int coarseTotal;
  int coarseOffset;
  int *coarseOffsets;
  int *coarseCounts;
  dfloat *invCoarseA;
  dfloat *xCoarse, *rhsCoarse;

  bool nullSpace;
  dfloat nullSpacePenalty;

  occa::device device;
  occa::stream defaultStream;
  occa::stream dataStream;  

  occa::memory o_x;
  occa::memory o_Ax;

  dfloat *rho;
  occa::memory o_rho;

  occa::kernel ellAXPYKernel;
  occa::kernel ellZeqAXPYKernel;
  occa::kernel ellJacobiKernel;
  occa::kernel cooAXKernel;
  occa::kernel scaleVectorKernel;
  occa::kernel vectorAddKernel;
  occa::kernel vectorAddKernel2;
  occa::kernel setVectorKernel;
  occa::kernel sumVectorKernel;
  occa::kernel addScalarKernel;
  occa::kernel dotStarKernel;
  occa::kernel simpleDotStarKernel;
  occa::kernel haloExtract;
  occa::kernel agg_interpolateKernel;
  occa::kernel innerProdKernel;
  occa::kernel vectorAddInnerProdKernel;
  occa::kernel kcycleCombinedOp1Kernel;
  occa::kernel kcycleCombinedOp2Kernel;
  occa::kernel vectorAddWeightedInnerProdKernel;
  occa::kernel kcycleWeightedCombinedOp1Kernel;
  occa::kernel kcycleWeightedCombinedOp2Kernel;

} parAlmond_t;

parAlmond_t *parAlmondInit(mesh_t *mesh, setupAide options);

void parAlmondAgmgSetup(parAlmond_t* parAlmond,
                       hlong* rowStarts,
                       dlong nnz,
                       hlong* Ai,
                       hlong* Aj,
                       dfloat* Avals,
                       bool nullSpace,
                       dfloat nullSpacePenalty);

void parAlmondPrecon(parAlmond_t* parAlmond, occa::memory o_x, occa::memory o_rhs);

int parAlmondFree(void* A);

#endif
