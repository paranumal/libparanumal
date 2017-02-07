#ifndef MESH_H 
#define MESH_H 1

#include <math.h>
#include <stdlib.h>
#include <occa.hpp>

#if 0
#define iint int
#define dfloat float
#define MPI_IINT MPI_INT
#define MPI_DFLOAT MPI_FLOAT
#define iintFormat "%d"
#define dfloatFormat "%f"
#define dfloatString "float"
#else
#define iint int
#define dfloat double
#define MPI_IINT MPI_INT
#define MPI_DFLOAT MPI_DOUBLE
#define iintFormat "%d"
#define dfloatFormat "%lf"
#define dfloatString "double"
#endif

#include "ogs_t.h"

typedef struct {

  iint dim;
  iint Nverts, Nfaces, NfaceVertices;

  iint Nnodes;
  dfloat *EX; // coordinates of vertices for each element
  dfloat *EY;
  dfloat *EZ;
  
  iint Nelements;
  iint *EToV; // element-to-vertex connectivity
  iint *EToE; // element-to-element connectivity
  iint *EToF; // element-to-(local)face connectivity
  iint *EToP; // element-to-partition/process connectivity
  iint *EToB; // element-to-boundary condition type

  // boundary faces
  iint NboundaryFaces; // number of boundary faces
  iint *boundaryInfo; // list of boundary faces (type, vertex-1, vertex-2, vertex-3)
  
  // MPI halo exchange info
  iint  totalHaloPairs;  // number of elements to be sent in halo exchange
  iint *haloElementList; // sorted list of elements to be sent in halo exchange
  iint *NhaloPairs;      // number of elements worth of data to send/recv
  iint  NhaloMessages;   // number of messages to send 
  
  void *haloSendRequests;
  void *haloRecvRequests;

  iint NinternalElements; // number of elements that can update without halo exchange
  iint NnotInternalElements; // number of elements that cannot update without halo exchange
  
  // NBN: streams / command queues
  occa::stream stream0, stream1;  
  
  // volumeGeometricFactors;
  dfloat *vgeo;
  iint Nvgeo;

  // second order volume geometric factors
  dfloat *ggeo;
  iint Nggeo;

  // volume node info 
  iint N, Np;
  dfloat *r, *s, *t;    // coordinates of local nodes
  dfloat *Dr, *Ds, *Dt; // collocation differentiation matrices
  dfloat *x, *y, *z;    // coordinates of physical nodes

  // indices of vertex nodes
  iint *vertexNodes;
  
  // quad specific quantity
  iint Nq, NqP;
  
  dfloat *D; // 1D differentiation matrix (for tensor-product)
  dfloat *gllz; // 1D GLL quadrature nodes
  dfloat *gllw; // 1D GLL quadrature weights

  // transform to/from eigenmodes of 1D laplacian (with built in weighting)
  dfloat *oasForward;
  dfloat *oasBack;
  dfloat *oasDiagOp;
  
  // transform to/from eigenmode of IPDG 1D laplacian
  dfloat *oasForwardDg;
  dfloat *oasBackDg;
  dfloat *oasDiagOpDg;
  
  // face node info
  iint Nfp;        // number of nodes per face
  iint *faceNodes; // list of element reference interpolation nodes on element faces
  iint *vmapM;     // list of volume nodes that are face nodes
  iint *vmapP;     // list of volume nodes that are paired with face nodes
  iint *mapP;     // list of surface nodes that are paired with -ve surface  nodes
  iint *faceVertices; // list of mesh vertices on each face

  dfloat *LIFT; // lift matrix

  iint   Nsgeo;
  dfloat *sgeo;

  // field info for PDE solver
  iint Nfields;
  dfloat *q;    // solution data array
  dfloat *rhsq; // right hand side data array
  dfloat *resq; // residual data array (for LSERK time-stepping)
  
  dfloat Lambda2; // square of penalty paramater used in constructing q^*

  // cubature
  iint cubNp;
  dfloat *cubr, *cubs, *cubt;    // coordinates of local nodes
  dfloat *cubx, *cuby, *cubz;    // coordinates of physical nodes
  dfloat *cubInterp; // interpolate from W&B to cubature nodes
  dfloat *cubProject; // projection matrix from cubature nodes to W&B nodes
  dfloat *cubDrW;    // 'r' weak differentiation matrix
  dfloat *cubDsW;    // 's' weak differentiation matrix
  dfloat *cubDtW;    // 't' weak differentiation matrix

  // c2 at cubature points (for wadg)
  dfloat *c2;

  // surface integration node info
  iint    intNfp;    // number of integration nodes on each face
  dfloat *intInterp; // interp from surface node to integration nodes
  dfloat *intLIFT;   // lift from surface integration nodes to W&B volume nodes
  dfloat *intx, *inty; // coordinates of suface integration nodes

  // Bernstein-Bezier info
  dfloat *VB, *invVB; // Bernstein Vandermonde matrices
  iint *D1ids, *D2ids, *D3ids; // Bernstein deriv matrix indices
  dfloat *Dvals; // Bernstein deriv matrix values
  dfloat *VBq, *PBq; // cubature interpolation/projection matrices
  dfloat *L0vals; // L0 values (L0 tridiagonal in 2D)
  iint *ELids; // lift reduction matrix indices
  dfloat *ELvals; // lift reduction matrix values
  iint max_EL_nnz; // max number of non-zeros per row of EL

  

  
  // time stepping info
  dfloat dt; // time step
  dfloat finalTime; // final time to run acoustics to
  iint   NtimeSteps;// number of time steps 
  iint   errorStep; // number of steps between error calculations
  iint   Nrk;
  dfloat rka[5], rkb[5], rkc[6];

  // ploting info for generating field vtu
  iint    plotNverts;    // number of vertices for each plot element
  iint    plotNp;        // number of plot nodes per element
  iint    plotNelements; // number of "plot elements" per element
  iint   *plotEToV;      // triangulation of plot nodes
  dfloat *plotR, *plotS, *plotT; // coordinates of plot nodes in reference element
  dfloat *plotInterp;    // warp & blend to plot node interpolation matrix

  // Boltzmann specific stuff
  dfloat RT, sqrtRT, tauInv; // need to remove this to ceedling

    // pml stuff
  iint    pmlNfields;
  //  iint    pmlNelements; // deprecated
  iint   *pmlElementList; // deprecated
  dfloat *pmlSigma;
  dfloat *pmlSigmaX;
  dfloat *pmlSigmaY;
  dfloat *pmlq;
  dfloat *pmlrhsq;
  dfloat *pmlresq;

  dfloat *pmlqx;    // x-pml data array
  dfloat *rhspmlqx; // right hand side data array
  dfloat *respmlqx; // residual data array (for LSERK time-stepping)
  dfloat *sigmax;

  dfloat *pmlqy;    // y-pml data array
  dfloat *rhspmlqy; // right hand side data array
  dfloat *respmlqy; // residual data array (for LSERK time-stepping)
  dfloat *sigmay;
  
  dfloat *pmlNT;    // time integrated relaxtion term
  dfloat *rhspmlNT; //
  dfloat *respmlNT; //

  
  
  // occa stuff
  occa::device device;
  occa::memory o_q, o_rhsq, o_resq;

  occa::memory o_Dr, o_Ds, o_Dt, o_LIFT;
  occa::memory o_DrT, o_DsT, o_DtT, o_LIFTT;

  occa::memory o_D; // tensor product differentiation matrix (for Hexes)
  
  occa::memory o_vgeo, o_sgeo;
  occa::memory o_vmapM, o_vmapP;
  
  occa::memory o_EToB, o_x, o_y, o_z;

  // cubature (for wadg)
  occa::memory o_intLIFTT, o_intInterpT, o_intx, o_inty;
  occa::memory o_cubDrWT, o_cubDsWT;
  occa::memory o_cubInterpT, o_cubProjectT;
  occa::memory o_invMc; // for comparison: inverses of weighted mass matrices
  occa::memory o_c2;

  // DG halo exchange info
  occa::memory o_haloElementList;
  occa::memory o_haloBuffer;

  occa::memory o_internalElementIds;
  occa::memory o_notInternalElementIds;
  
  // Bernstein-Bezier occa arrays
  occa::memory o_D1ids, o_D2ids, o_D3ids, o_Dvals; // Bernstein deriv matrix indices
  occa::memory o_VBq, o_PBq; // cubature interpolation/projection matrices
  occa::memory o_L0vals, o_ELids, o_ELvals; 


  // pml vars
  occa::memory o_sigmax, o_sigmay;

  iint pmlNelements;
  iint nonPmlNelements;
  occa::memory o_pmlElementIds;
  occa::memory o_nonPmlElementIds;

  occa::memory o_pmlqx, o_rhspmlqx, o_respmlqx;
  occa::memory o_pmlqy, o_rhspmlqy, o_respmlqy;
  occa::memory o_pmlNT, o_rhspmlNT, o_respmlNT;
  
  occa::memory o_pmlElementList;
  occa::memory o_pmlSigmaX, o_pmlSigmaY;
  
  occa::memory o_pmlrhsq, o_pmlresq, o_pmlq;


  
  // CG gather-scatter info
  void *gsh; // gslib struct pointer

  iint *gatherLocalIds; // local index of rank/gather id sorted list of nodes
  iint *gatherBaseIds;  // gather index of ""
  iint *gatherBaseRanks; // base rank
  iint *gatherMaxRanks;  // maximum rank connected to each sorted node
  iint *gatherHaloFlags;  // maximum rank connected to each sorted node

  iint NuniqueBases; // number of unique bases on this rank
  occa::memory o_gatherNodeOffsets; // list of offsets into gatherLocalNodes for start of base
  occa::memory o_gatherLocalNodes; // indices of local nodes collected by base node
  occa::memory o_gatherTmp; // temporary array to store base values gathered locally

  iint NnodeHalo; // number of halo bases on this rank
  occa::memory o_nodeHaloIds;  // indices of halo base nodes after initial local gather
  occa::memory o_subGatherTmp; // temporary DEVICE array to store halo base values prior to DEVICE>HOST copy
  dfloat        *subGatherTmp; // temporary HALO array

  occa::memory o_ggeo; // second order geometric factors
  occa::memory o_projectL2; // local weights for projection.

  
  occa::kernel volumeKernel;
  occa::kernel surfaceKernel;
  occa::kernel updateKernel;
  occa::kernel haloExtractKernel;
  occa::kernel partialSurfaceKernel;
  
  occa::kernel gatherKernel;
  occa::kernel scatterKernel;

  occa::kernel getKernel;
  occa::kernel putKernel;

  occa::kernel AxKernel;
  occa::kernel innerProductKernel;
  occa::kernel weightedInnerProduct1Kernel;
  occa::kernel weightedInnerProduct2Kernel;
  occa::kernel scaledAddKernel;
  occa::kernel dotMultiplyKernel;
  occa::kernel dotDivideKernel;

  occa::kernel gradientKernel;
  occa::kernel ipdgKernel;

  occa::kernel relaxationKernel;
  
  occa::kernel pmlKernel; // deprecated
  occa::kernel pmlVolumeKernel;
  occa::kernel pmlSurfaceKernel;
  occa::kernel pmlUpdateKernel;

  
}mesh_t;

// serial sort
void mysort(iint *data, iint N, const char *order);

// sort entries in an array in parallel
void parallelSort(iint N, void *vv, size_t sz,
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

#endif

