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

#ifndef MESH_H
#define MESH_H 1

#include <unistd.h>

#include "mpi.h"
#include <math.h>
#include <stdlib.h>
#include <occa.hpp>

#include "types.h"
#include "ogs_t.h"

#include "timer.h"

#include "setupAide.hpp"

#define TRIANGLES 3
#define QUADRILATERALS 4
#define TETRAHEDRA 6
#define HEXAHEDRA 12

typedef struct {

  MPI_Comm comm;
  int rank, size; // MPI rank and size (process count)
  
  int dim;
  int Nverts, Nfaces, NfaceVertices;

  hlong Nnodes;
  dfloat *EX; // coordinates of vertices for each element
  dfloat *EY;
  dfloat *EZ;

  dlong Nelements;
  hlong *EToV; // element-to-vertex connectivity
  dlong *EToE; // element-to-element connectivity
  int   *EToF; // element-to-(local)face connectivity
  int   *EToP; // element-to-partition/process connectivity
  int   *EToB; // element-to-boundary condition type

  int *elementInfo; //type of element

  // boundary faces
  hlong NboundaryFaces; // number of boundary faces
  hlong *boundaryInfo; // list of boundary faces (type, vertex-1, vertex-2, vertex-3)

  // MPI halo exchange info
  dlong  totalHaloPairs;  // number of elements to be sent in halo exchange
  dlong *haloElementList; // sorted list of elements to be sent in halo exchange
  int *NhaloPairs;      // number of elements worth of data to send/recv
  int  NhaloMessages;     // number of messages to send

  void *haloSendRequests;
  void *haloRecvRequests;

  dlong NinternalElements; // number of elements that can update without halo exchange
  dlong NnotInternalElements; // number of elements that cannot update without halo exchange

  //list of fair pairs
  dlong NfacePairs;
  dlong *EToFPairs;
  dlong *FPairsToE;
  int *FPairsToF;

  // NBN: streams / command queues
  occa::stream stream0, stream1;

  // volumeGeometricFactors;
  dlong Nvgeo;
  dfloat *vgeo;

  // second order volume geometric factors
  dlong Nggeo;
  dfloat *ggeo;

  // volume node info
  int N, Np;
  dfloat *r, *s, *t;    // coordinates of local nodes
  dfloat *Dr, *Ds, *Dt; // collocation differentiation matrices
  dfloat *Dmatrices;
  dfloat *MM, *invMM;           // reference mass matrix
  dfloat *Srr,*Srs, *Srt; //element stiffness matrices
  dfloat *Ssr,*Sss, *Sst;
  dfloat *Str,*Sts, *Stt;
  dfloat *Smatrices;
  int maxNnzPerRow;
  dfloat *x, *y, *z;    // coordinates of physical nodes
  
  dfloat sphereRadius;  // for Quad3D 
  

  // indices of vertex nodes
  int *vertexNodes;

  // quad specific quantity
  int Nq, NqP, NpP;

  dfloat *D; // 1D differentiation matrix (for tensor-product)
  dfloat *gllz; // 1D GLL quadrature nodes
  dfloat *gllw; // 1D GLL quadrature weights

  int gjNq;
  dfloat *gjr,*gjw; // 1D nodes and weights for Gauss Jacobi quadature
  dfloat *gjI,*gjD; // 1D GLL to Gauss node interpolation and differentiation matrices
  dfloat *gjD2;     // 1D GJ to GJ node differentiation

  // transform to/from eigenmodes of 1D laplacian (with built in weighting)
  dfloat *oasForward;
  dfloat *oasBack;
  dfloat *oasDiagOp;

  // transform to/from eigenmode of IPDG 1D laplacian
  dfloat *oasForwardDg;
  dfloat *oasBackDg;
  dfloat *oasDiagOpDg;

  //rotated node ids
  int *rmapP;

  //reference patch inverse (for OAS precon)
  dfloat *invAP;

  // face node info
  int Nfp;        // number of nodes per face
  int *faceNodes; // list of element reference interpolation nodes on element faces
  dlong *vmapM;     // list of volume nodes that are face nodes
  dlong *vmapP;     // list of volume nodes that are paired with face nodes
  dlong *mapP;     // list of surface nodes that are paired with -ve surface  nodes
  int *faceVertices; // list of mesh vertices on each face

  dfloat *LIFT; // lift matrix
  dfloat *FMM;  // Face Mass Matrix
  dfloat *sMT; // surface mass (MM*LIFT)^T

  dlong   Nsgeo;
  dfloat *sgeo;

  // field info for PDE solver
  int Nfields;
  dfloat *q;    // solution data array
  dfloat *fQM, *fQP; //solution trace arrays
  dfloat *rhsq, *rhsq2, *rhsq3; // right hand side data array
  dfloat *resq; // residual data array (for LSERK time-stepping)

  dfloat Lambda2; // square of penalty paramater used in constructing q^*

  // cubature
  int cubNp, cubNfp, cubNq;
  dfloat *cubr, *cubs, *cubt, *cubw; // coordinates and weights of local cubature nodes
  dfloat *cubx, *cuby, *cubz;    // coordinates of physical nodes
  dfloat *cubInterp; // interpolate from W&B to cubature nodes
  dfloat *cubProject; // projection matrix from cubature nodes to W&B nodes
  dfloat *cubDW;     // 1D weak differentiation matrix
  dfloat *cubDrW;    // 'r' weak differentiation matrix
  dfloat *cubDsW;    // 's' weak differentiation matrix
  dfloat *cubDtW;    // 't' weak differentiation matrix
  dfloat *cubDWmatrices;

  dfloat *cubvgeo;  //volume geometric data at cubature points
  dfloat *cubsgeo;  //surface geometric data at cubature points

  // c2 at cubature points (for wadg)
  dfloat *c2;

  //source injection
  dfloat *sourceq;
  dfloat sourceX0, sourceY0, sourceZ0, sourceT0, sourceC2, sourceFreq;
  int sourceNelements;
  dlong *MRABsourceNelements;
  dlong *sourceElements;

  // surface integration node info
  int    intNfp;    // number of integration nodes on each face
  dfloat *intInterp; // interp from surface node to integration nodes
  dfloat *intLIFT;   // lift from surface integration nodes to W&B volume nodes
  dfloat *intx, *inty, *intz; // coordinates of suface integration nodes

  // Bernstein-Bezier info
  dfloat *VB, *invVB; // Bernstein Vandermonde matrices
  dfloat *BBMM;
  dfloat *invVB1D, *invVB2D;
  int *D0ids, *D1ids, *D2ids, *D3ids; // Bernstein deriv matrix indices
  dfloat *Dvals; // Bernstein deriv matrix values
  int *D0Tids, *D1Tids, *D2Tids, *D3Tids; // Bernstein transpose deriv matrix indices
  dfloat *DTvals; // Bernstein transpose deriv matrix values
  dfloat *VBq, *PBq; // cubature interpolation/projection matrices
  int *L0ids; // L0 matrix ids
  dfloat *L0vals; // L0 values (L0 tridiagonal in 2D)
  int *ELids; // lift reduction matrix indices
  dfloat *ELvals; // lift reduction matrix values
  int max_EL_nnz; // max number of non-zeros per row of EL
  int *BBRaiseids; //Bernstein elevate matrix indices
  dfloat *BBRaiseVals; //Bernstein elevate matrix values
  dfloat *BBLower; //Berstein projection matrix.

  //degree raising and lowering interpolation matrices
  dfloat *interpRaise;
  dfloat *interpLower;

  //sparse basis info
  dfloat *sparseV, *invSparseV;
  dfloat *sparseMM;
  int* FaceModes;
  int SparseNnzPerRow;
  int SparseNnzPerRowNonPadded;
  int *sparseStackedNZ;
  dfloat *sparseSrrT;
  dfloat *sparseSrsT;
  dfloat *sparseSssT;
  int *Ind;

  dlong *mmapM, *mmapP; 
  int   *mmapS;
  dfloat *mapSgn;

  // time stepping info
  dfloat dt; // time step
  dfloat startTime ; // Start Time
  dfloat finalTime; // final time to run acoustics to
  int   NtimeSteps;// number of time steps
  int   errorStep; // number of steps between error calculations
  int   Nrk;
  dfloat rka[5], rkb[5], rkc[6]; // AK: deprecated

  // MRAB,SAAB coefficients
  dfloat mrab[3], mrabb[3], saab[3], saabexp; // AK: deprecated 
  int MRABNlevels;
  int *MRABlevel;
  dlong *MRABNelements, *MRABNhaloElements;
  dlong **MRABelementIds, **MRABhaloIds;
  int *MRABshiftIndex;

  dlong *MRABpmlNelements, *MRABpmlNhaloElements;
  dlong **MRABpmlElementIds, **MRABpmlIds;
  dlong **MRABpmlHaloElementIds, **MRABpmlHaloIds;

  dlong pmlNelements, nonPmlNelements;
  dlong *nonPmlElementIds, *pmlElementIds, *pmlIds;  
  int shiftIndex;

  dfloat dtfactor; //Delete later for script run 
  dfloat maxErrorBoltzmann;

  //LSIMEX-BOLTZMANN coefficients, simplified for efficient implementation
  dfloat LSIMEX_B[4], LSIMEX_C[4], LSIMEX_ABi[4], LSIMEX_ABe[4], LSIMEX_Ad[4];
  dfloat *MRSAAB_A, *MRSAAB_B, *MRSAAB_C, *MRAB_A, *MRAB_B, *MRAB_C;
  dfloat RK_A[5][5], RK_B[5], RK_C[5], SARK_A[5][5], SARK_B[5], SARK_C[5]; 
  int Nimex, Nrhs; // deprecated
  
  dfloat *rkq, *rkrhsq, *rkerr; // deprecated, AK. 
  dfloat *errtmp;
  dfloat rkC[7], rkA[7*7], rkE[7];

  occa::memory o_rkq, o_rkrhsq, o_rkerr; // deprecated, AK.
  occa::memory o_errtmp;
  occa::memory o_rkA, o_rkE;

  // ploting info for generating field vtu
  int    plotNverts;    // number of vertices for each plot element
  int    plotNp;        // number of plot nodes per element
  int    plotNelements; // number of "plot elements" per element
  int   *plotEToV;      // triangulation of plot nodes
  dfloat *plotR, *plotS, *plotT; // coordinates of plot nodes in reference element
  dfloat *plotInterp;    // warp & blend to plot node interpolation matrix

  int *contourEToV;
  dfloat *contourVX, *contourVY, *contourVZ;
  dfloat *contourInterp, *contourInterp1, *contourFilter; 

  //SEMFEM data
  int NpFEM, NelFEM;
  int *FEMEToV;
  dfloat *rFEM, *sFEM, *tFEM;
  dfloat *SEMFEMInterp;

  occa::memory o_SEMFEMInterp;
  occa::memory o_SEMFEMAnterp;

  // Boltzmann specific stuff
  dfloat RT, sqrtRT, tauInv, Ma, Re; // Deprecated: AK

  // pml stuff
  int    pmlNfields;
  //  dlong    pmlNelements; // deprecated
  dlong   *pmlElementList; // deprecated

  int Ntscale; // Will be removed, for time accuracy test
  
  dfloat *pmlBetaX, *pmlBetaY; // deprecated
  dfloat *pmlSigma;
  dfloat *pmlSigmaX;
  dfloat *pmlSigmaY;
  dfloat *pmlSigmaZ;

  dfloat *pmlq;
  dfloat *pmlqx;
  dfloat *pmlqy;
  dfloat *pmlqz;

  dfloat *pmlrhsq;
  dfloat *pmlrhsqx;
  dfloat *pmlrhsqy;
  dfloat *pmlrhsqz;

  dfloat *pmlresq;
  dfloat *pmlresqx;
  dfloat *pmlresqy;
  dfloat *pmlresqz;



  dfloat *invTau; // deprecated in Boltzmann


  // AK: Remove the below definition after fixing MRAB, only single rate uses 
  // dfloat *pmlqx;    // x-pml data array
  dfloat *rhspmlqx; // right hand side data array
  dfloat *respmlqx; // residual data array (for LSERK time-stepping)
  dfloat *sigmax;

  // dfloat *pmlqy;    // y-pml data array
  dfloat *rhspmlqy; // right hand side data array
  dfloat *respmlqy; // residual data array (for LSERK time-stepping)
  dfloat *sigmay;

  // dfloat *pmlqz;    // Z-pml data array
  dfloat *rhspmlqz; // right hand side data array
  dfloat *respmlqz; // residual data array (for LSERK time-stepping)
  dfloat *sigmaz;

  //dfloat *pmlq;    // Z-pml data array
  dfloat *rhspmlq; // right hand side data array
  dfloat *respmlq; // residual data array (for LSERK time-stepping)

  // Probe Data
  int probeN, probeNTotal; 
  dfloat *probeR, *probeS, *probeT;
  // dfloat *probeX, *probeY, *probeZ;  
  dlong *probeElementIds, *probeIds;  
  dfloat *probeI; 





  // occa stuff
  occa::device device;

  occa::stream defaultStream;
  occa::stream dataStream;

  occa::memory o_q, o_rhsq, o_resq, o_fQM, o_fQP;

  occa::memory o_Dr, o_Ds, o_Dt, o_LIFT, o_MM;
  occa::memory o_DrT, o_DsT, o_DtT, o_LIFTT;
  occa::memory o_Dmatrices;
  occa::memory o_FMMT;
  occa::memory o_sMT;

  occa::memory o_D; // tensor product differentiation matrix (for Hexes)
  occa::memory o_SrrT, o_SrsT, o_SrtT; //element stiffness matrices
  occa::memory o_SsrT, o_SssT, o_SstT;
  occa::memory o_Srr, o_Srs, o_Srt, o_Sss, o_Sst, o_Stt; // for char4-based kernels
  occa::memory o_Smatrices;
  occa::memory o_IndT, o_IndTchar;
  occa::memory o_India, o_Indja;
  occa::memory o_StrT, o_StsT, o_SttT;
  occa::memory o_Ind; // for sparse index storage

  occa::memory o_vgeo, o_sgeo;
  occa::memory o_vmapM, o_vmapP, o_mapP;

  occa::memory o_rmapP;

  occa::memory o_EToE, o_EToF, o_EToB, o_x, o_y, o_z;

  occa::memory o_EToFPairs, o_FPairsToE, o_FPairsToF;

  // cubature (for wadg)
  occa::memory o_intLIFTT, o_intInterpT, o_intx, o_inty, o_intz;
  occa::memory o_cubDWT;
  occa::memory o_cubDrWT, o_cubDsWT, o_cubDtWT;
  occa::memory o_cubDWmatrices;
  occa::memory o_cubInterpT, o_cubProjectT;
  occa::memory o_invMc; // for comparison: inverses of weighted mass matrices

  occa::memory o_cubvgeo, o_cubsgeo;

  occa::memory o_c2;

  //MRAB element lists
  occa::memory *o_MRABelementIds;
  occa::memory *o_MRABhaloIds;
  occa::memory *o_MRABpmlElementIds;
  occa::memory *o_MRABpmlIds;
  occa::memory *o_MRABpmlHaloElementIds;
  occa::memory *o_MRABpmlHaloIds;


  // DG halo exchange info
  occa::memory o_haloElementList;
  occa::memory o_haloBuffer;

  occa::memory o_internalElementIds;
  occa::memory o_notInternalElementIds;

  // Bernstein-Bezier occa arrays
  occa::memory o_BBMM;
  occa::memory o_D0ids, o_D1ids, o_D2ids, o_D3ids, o_Dvals; // Bernstein deriv matrix indices
  occa::memory o_packedDids; // char4 packed increments (D1ids-D0ids)

  occa::memory o_invVB1DT, o_invVB2DT;
  occa::memory o_VBq, o_PBq; // cubature interpolation/projection matrices
  occa::memory o_L0ids, o_L0vals, o_ELids, o_ELvals;

  /* sparse basis occa arrays */
  occa::memory o_sparseStackedNZ;
  occa::memory o_sparseSrrT;
  occa::memory o_sparseSrsT;
  occa::memory o_sparseSssT;
  occa::memory o_mapSgn;

  // pml vars
  occa::memory o_sigmax, o_sigmay, o_sigmaz; // AK: deprecated


  occa::memory o_pmlElementIds;
  occa::memory o_nonPmlElementIds;
  occa::memory o_pmlIds;

  occa::memory o_pmlElementList;
  
  // Deprecated in Boltzmann
  occa::memory o_pmlSigmaX, o_pmlSigmaY, o_pmlSigmaZ;
  occa::memory o_pmlBetaX, o_pmlBetaY; // deprecated
  occa::memory o_pmlq, o_pmlrhsq, o_pmlresq ; 
  occa::memory o_pmlqx,o_pmlqy, o_pmlqz; 
  occa::memory o_pmlrhsqx, o_pmlrhsqy, p_pmlrhsqz;
  occa::memory o_pmlresqx, o_pmlresqy, p_pmlresqz;


  // occa::memory o_rhspmlqx, o_respmlqx; 
  // occa::memory o_rhspmlqy, o_respmlqy;
  // occa::memory o_rhspmlqz, o_respmlqz;
  occa::memory o_pmlNT, o_rhspmlNT, o_respmlNT; // deprecated !

  // Boltzmann SARK extra storage for exponential update
  // occa::memory o_resqex;

  // Boltzmann SAAB 3th order storage: respmlqx, qy, nt and q not used
  occa::memory o_expsigmax, o_expsigmay; // deprecated
  occa::memory o_rhsq2,     o_rhsq3;     // deprecated
  occa::memory o_rhspmlqx2, o_rhspmlqx3; // deprecated
  occa::memory o_rhspmlqy2, o_rhspmlqy3; // deprecated
  occa::memory o_rhspmlNT2, o_rhspmlNT3; // deprecated
  // Deprecated
  occa::memory o_qY,   o_qZ,   o_qS;
  occa::memory o_qYx,  o_qZx,  o_qSx;
  occa::memory o_qYy,  o_qZy,  o_qSy;








  // AK: Remove this stuff, rename single rate files
  occa::memory o_rhspmlq,   o_respmlq; // 3D LSERK
  occa::memory o_pmlqold,  o_rhspmlq2,  o_rhspmlq3; // 3D Semianalytic
  occa::memory o_pmlqY, o_pmlqS; // 3D IMEX




  // CG gather-scatter info
  void *gsh, *hostGsh; // gslib struct pointer
  ogs_t *ogs; //occa gs pointer

  hlong *globalIds;
  int *globalOwners;
  int *globalHaloFlags;

  dlong *gatherLocalIds; // local index of rank/gather id sorted list of nodes
  hlong *gatherBaseIds;  // gather index of ""
  int *gatherBaseRanks; // base rank
  dlong *gatherOffsets; 
  int *gatherMaxRanks;  // maximum rank connected to each sorted node
  int *gatherHaloFlags;  // maximum rank connected to each sorted node
  hlong *gatherGlobalStarts;

  dlong NuniqueBases; // number of unique bases on this rank
  occa::memory o_gatherNodeOffsets; // list of offsets into gatherLocalNodes for start of base
  occa::memory o_gatherLocalNodes; // indices of local nodes collected by base node
  occa::memory o_gatherTmp; // temporary array to store base values gathered locally

  dlong NnodeHalo; // number of halo bases on this rank
  occa::memory o_nodeHaloIds;  // indices of halo base nodes after initial local gather
  occa::memory o_subGatherTmp; // temporary DEVICE array to store halo base values prior to DEVICE>HOST copy
  dfloat        *subGatherTmp; // temporary HALO array

  occa::memory o_ggeo; // second order geometric factors
  occa::memory o_projectL2; // local weights for projection.

  occa::kernel volumeKernel;
  occa::kernel surfaceKernel;
  occa::kernel updateKernel;
  occa::kernel traceUpdateKernel;
  occa::kernel haloExtractKernel;
  occa::kernel partialSurfaceKernel;
  

  // Just for test will be deleted after temporal testsAK
  occa::kernel RKupdateKernel;
  occa::kernel RKpmlUpdateKernel;


  occa::kernel gatherKernel;
  occa::kernel scatterKernel;
  occa::kernel gatherScatterKernel;

  occa::kernel getKernel;
  occa::kernel putKernel;

  occa::kernel sumKernel;
  occa::kernel addScalarKernel;

  occa::kernel AxKernel;
  occa::kernel innerProductKernel;
  occa::kernel weightedInnerProduct1Kernel;
  occa::kernel weightedInnerProduct2Kernel;
  occa::kernel scaledAddKernel;
  occa::kernel dotMultiplyKernel;
  occa::kernel dotDivideKernel;

  occa::kernel gradientKernel;
  occa::kernel ipdgKernel;

  occa::kernel maskKernel;

  // Boltzmann Specific Kernels
  occa::kernel relaxationKernel;
  occa::kernel pmlRelaxationKernel;
  // //Boltzmann Imex Kernels
  

  //Deprecated
  occa::kernel implicitUpdateKernel;
  occa::kernel pmlImplicitUpdateKernel;
  occa::kernel implicitSolveKernel;
  occa::kernel pmlImplicitSolveKernel;
  //
  occa::kernel residualUpdateKernel;
  occa::kernel pmlResidualUpdateKernel;



  occa::kernel pmlKernel;
  occa::kernel pmlVolumeKernel;
  occa::kernel pmlSurfaceKernel;
  occa::kernel pmlUpdateKernel;
  occa::kernel pmlTraceUpdateKernel;


  // Experimental Time Steppings for Boltzmann
#if 1
  occa::kernel updateStageKernel;
  occa::kernel pmlUpdateStageKernel;
  //occa::kernel updateStageKernel33;
  //
  occa::memory o_rhsq4, o_rhsq5, o_invTau;
  occa::memory o_qold , o_pmlqxold, o_pmlqyold,o_pmlNTold;
  // SARK extra coefficients for Boltzmann Solver
  dfloat sarka[5][5], sarkb[5], sarke[6], sarkra[5], sarkrb[5]; // exponential update terms, better to hold
  dfloat sarkpmla[5][5], sarkpmlb[5], sarkpmle[6];
  dfloat rk3a[3][3], rk3b[3], rk3c[3];
  dfloat rk4a[5][5], rk4b[5];
  dfloat lserk3a[3], lserk3b[3], lserk3c[4];

#endif



}mesh_t;

// serial sort
void mysort(hlong *data, int N, const char *order);

// sort entries in an array in parallel
void parallelSort(int size, int rank, MPI_Comm comm,
		  int N, void *vv, size_t sz,
		  int (*compare)(const void *, const void *),
		  void (*match)(void *, void *)
		  );

#define mymax(a,b) (((a)>(b))?(a):(b))
#define mymin(a,b) (((a)<(b))?(a):(b))

/* hash function */
unsigned int hash(const unsigned int value) ;

/* dimension independent mesh operations */
void meshConnect(mesh_t *mesh);

/* build parallel face connectivity */
void meshParallelConnect(mesh_t *mesh);

/* build global connectivity in parallel */
void meshParallelConnectNodes(mesh_t *mesh);

void meshHaloSetup(mesh_t *mesh);

/* extract whole elements for the halo exchange */
void meshHaloExtract(mesh_t *mesh, size_t Nbytes, void *sourceBuffer, void *haloBuffer);

void meshHaloExchange(mesh_t *mesh,
    size_t Nbytes,         // message size per element
    void *sourceBuffer,
    void *sendBuffer,    // temporary buffer
    void *recvBuffer);

void meshHaloExchangeStart(mesh_t *mesh,
    size_t Nbytes,       // message size per element
    void *sendBuffer,    // temporary buffer
    void *recvBuffer);


void meshHaloExchangeFinish(mesh_t *mesh);

// print out parallel partition i
void meshPartitionStatistics(mesh_t *mesh);

// build element-boundary connectivity
void meshConnectBoundary(mesh_t *mesh);

// squeeze gaps out of a globalNumbering of local nodes (arranged in NpNum blocks
void meshParallelConsecutiveGlobalNumbering(mesh_t *mesh,
                                            dlong Nnum,
                                            hlong *globalNumbering, 
                                            int *globalOwners, 
                                            hlong *globalStarts);

ogs_t *meshParallelGatherScatterSetup(mesh_t *mesh,
                                      dlong Nlocal,
                                      dlong *gatherLocalIds,
                                      hlong *gatherBaseIds,
                                      int *gatherBaseRanks,
                                      int  *gatherHaloFlags,
                                      int verbose);

void meshParallelGatherScatter(mesh_t *mesh, ogs_t *ogs, occa::memory &o_v);
void meshParallelGather(mesh_t *mesh, ogs_t *ogs, occa::memory &o_v, occa::memory &o_gv);
void meshParallelScatter(mesh_t *mesh, ogs_t *ogs, occa::memory &o_v, occa::memory &o_sv);


void occaTimerTic(occa::device device,std::string name);
void occaTimerToc(occa::device device,std::string name);

extern "C"
{
  void *gsParallelGatherScatterSetup(MPI_Comm comm, dlong Ngather, hlong *gatherIds, int verbose);
  void gsParallelGatherScatter(void *gsh, void *v, const char *type, const char *op);
  void gsVecParallelGatherScatter(void *gsh, void *v, int k, const char *type, const char *op);
  void gsParallelGatherScatterDestroy(void *gsh);

  void * xxtSetup(uint num_local_rows,
      void* row_ids,
      uint nnz,
      void*   A_i,
      void*   A_j,
      void* A_vals,
      int null_space,
      const char* inttype,
      const char* floattype);

  void xxtSolve(void* x,
      void* A,
      void* rhs);

  void xxtFree(void* A) ;
}

extern "C"
{
  void dgesv_ ( int     *N, int     *NRHS, double  *A,
                int     *LDA,
                int     *IPIV, 
                double  *B,
                int     *LDB,
                int     *INFO );

  void sgesv_(int *N, int *NRHS,float  *A, int *LDA, int *IPIV, float  *B, int *LDB,int *INFO);

  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
  void dgeev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, double *WR, double *WI,
              double *VL, int *LDVL, double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO );
  
  double dlange_(char *NORM, int *M, int *N, double *A, int *LDA, double *WORK);
  void dgecon_(char *NORM, int *N, double *A, int *LDA, double *ANORM,
                double *RCOND, double *WORK, int *IWORK, int *INFO );
}

void readDfloatArray(FILE *fp, const char *label, dfloat **A, int *Nrows, int* Ncols);
void readIntArray   (FILE *fp, const char *label, int **A   , int *Nrows, int* Ncols);

void meshApplyElementMatrix(mesh_t *mesh, dfloat *A, dfloat *q, dfloat *Aq);

void matrixInverse(int N, dfloat *A);
dfloat matrixConditionNumber(int N, dfloat *A);

void occaDeviceConfig(mesh_t *mesh, setupAide &newOptions);

void *occaHostMallocPinned(occa::device &device, size_t size, void *source, occa::memory &mem);

#endif

