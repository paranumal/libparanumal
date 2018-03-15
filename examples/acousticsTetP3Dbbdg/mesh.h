#ifndef MESH_H 
#define MESH_H 1

#include <math.h>
#include <stdlib.h>
#include <occa.hpp>

#if 0
#define dfloat float
#define MPI_DFLOAT MPI_FLOAT
#define dfloatFormat "%f"
#define dfloatString "float"
#else
#define dfloat double
#define MPI_DFLOAT MPI_DOUBLE
#define dfloatFormat "%lf"
#define dfloatString "double"
#endif


typedef struct {

  int dim;
  int Nverts, Nfaces, NfaceVertices;

  int Nnodes;
  dfloat *EX; // coordinates of vertices for each element
  dfloat *EY;
  dfloat *EZ;
  
  int Nelements;
  int *EToV; // element-to-vertex connectivity
  int *EToE; // element-to-element connectivity
  int *EToF; // element-to-(local)face connectivity
  int *EToP; // element-to-partition/process connectivity
  int *EToB; // element-to-boundary condition type

  int *elementInfo; //type of element

  // boundary faces
  int NboundaryFaces; // number of boundary faces
  int *boundaryInfo; // list of boundary faces (type, vertex-1, vertex-2, vertex-3)
  
  // MPI halo exchange info
  int  totalHaloPairs;  // number of elements to be sent in halo exchange
  int *haloElementList; // sorted list of elements to be sent in halo exchange
  int *NhaloPairs;      // number of elements worth of data to send/recv
  int  NhaloMessages;   // number of messages to send 
  
  void *haloSendRequests;
  void *haloRecvRequests;

  int NinternalElements; // number of elements that can update without halo exchange
  int NnotInternalElements; // number of elements that cannot update without halo exchange
  
  // NBN: streams / command queues
  occa::stream stream0, stream1;  
  
  // volumeGeometricFactors;
  dfloat *vgeo;
  int Nvgeo;

  // second order volume geometric factors
  dfloat *ggeo;
  int Nggeo;

  // volume node info 
  int NMax, NpMax, NfpMax;
  int *N, *Np;
  int *NelOrder, **NelList;

  dfloat **r, **s, **t;    // coordinates of local nodes
  dfloat **Dr,**Ds, **Dt; // collocation differentiation matrices
  dfloat **MM;           // reference mass matrix
  dfloat *x, *y, *z;    // coordinates of physical nodes

  // indices of vertex nodes
  int *vertexNodes;
  
  // quad specific quantity
  int Nq, NqP;
  
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
  int *Nfp;        // number of nodes per face
  int **faceNodes; // list of element reference interpolation nodes on element faces
  int *vmapM;     // list of volume nodes that are face nodes
  int *vmapP;     // list of volume nodes that are paired with face nodes
  int *mapP;     // list of surface nodes that are paired with -ve surface  nodes
  int *faceVertices; // list of mesh vertices on each face

  dfloat **LIFT; // lift matrix

  int   Nsgeo;
  dfloat *sgeo;

  // field info for PDE solver
  int Nfields;
  dfloat *q;    // solution data array
  dfloat *fQM, *fQP; //solution trace array
  dfloat *rhsq; // right hand side data array
  dfloat *resq; // residual data array (for LSERK time-stepping)
  
  dfloat Lambda2; // square of penalty paramater used in constructing q^*

  // cubature
  int *cubNp, cubNpMax;
  dfloat **cubr, **cubs, **cubt;    // coordinates of local nodes
  dfloat **cubx, **cuby, **cubz;    // coordinates of physical nodes
  dfloat **cubInterp; // interpolate from W&B to cubature nodes
  dfloat **cubProject; // projection matrix from cubature nodes to W&B nodes
  dfloat **cubDrW;    // 'r' weak differentiation matrix
  dfloat **cubDsW;    // 's' weak differentiation matrix
  dfloat **cubDtW;    // 't' weak differentiation matrix

  // c2 at cubature points (for wadg)
  dfloat *c2;

  //source injection
  dfloat *sourceq;
  dfloat sourceX0, sourceY0, sourceZ0, sourceT0, sourceC2, sourceFreq;
  int sourceNelements, *MRABsourceNelements;
  int *sourceElements;

  // surface integration node info
  int    *intNfp;    // number of integration nodes on each face
  dfloat **intInterp; // interp from surface node to integration nodes
  dfloat **intLIFT;   // lift from surface integration nodes to W&B volume nodes
  dfloat *intx, *inty; // coordinates of suface integration nodes

  // Bernstein-Bezier info
  dfloat **VB, **invVB; // Bernstein Vandermonde matrices
  dfloat **invVB1D, **invVB2D;
  int **D0ids, **D1ids, **D2ids, **D3ids; // Bernstein deriv matrix indices
  dfloat **Dvals; // Bernstein deriv matrix values
  dfloat **VBq, **PBq; // cubature interpolation/projection matrices
  dfloat **VBplot;
  int **L0ids; // L0 matrix ids
  dfloat **L0vals; // L0 values (L0 tridiagonal in 2D)
  int **ELids; // lift reduction matrix indices
  dfloat **ELvals; // lift reduction matrix values
  int *max_EL_nnz; // max number of non-zeros per row of EL
  int **BBRaiseids; //Bernstein elevate matrix indices
  dfloat **BBRaiseVals; //Bernstein elevate matrix values
  dfloat **BBLower; //Berstein projection matrix.
  
  // time stepping info
  dfloat dt; // time step
  dfloat finalTime; // final time to run acoustics to
  int   NtimeSteps;// number of time steps 
  int   errorStep; // number of steps between error calculations
  int   Nrk;
  dfloat rka[5], rkb[5], rkc[6];

  // MRAB lists
  int MRABNlevels;
  int *MRABlevel;
  int *MRABNelements, *MRABNhaloElements;
  int **MRABelementIds, **MRABhaloIds;
  int **MRABNelP, **MRABNhaloEleP;
  int ***MRABelIdsP, ***MRABhaloIdsP;
  int *MRABshiftIndex;

  int *MRABpmlNelements, *MRABpmlNhaloElements;
  int **MRABpmlElementIds, **MRABpmlIds;
  int **MRABpmlHaloElementIds, **MRABpmlHaloIds;

  int **MRABpmlNelP, **MRABpmlNhaloEleP;
  int ***MRABpmlElIdsP, ***MRABpmlIdsP;
  int ***MRABpmlHaloEleIdsP, ***MRABpmlHaloIdsP;

  // ploting info for generating field vtu
  int    plotNverts;    // number of vertices for each plot element
  int    *plotNp;        // number of plot nodes per element
  int    *plotNelements; // number of "plot elements" per element
  int   **plotEToV;      // triangulation of plot nodes
  dfloat **plotR, **plotS, **plotT; // coordinates of plot nodes in reference element
  dfloat **plotInterp;    // warp & blend to plot node interpolation matrix

  // Boltzmann specific stuff
  dfloat RT, sqrtRT, tauInv; // need to remove this to ceedling

    // pml stuff
  int    pmlNfields;
  //  int    pmlNelements; // deprecated
  int   *pmlElementList; // deprecated
  dfloat *pmlSigma;
  dfloat *pmlSigmaX;
  dfloat *pmlSigmaY;
  dfloat *pmlSigmaZ;
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
  occa::memory o_fQM, o_fQP;

  occa::memory o_N;

  occa::memory *o_NelList;

  occa::memory o_Dr, o_Ds, o_Dt, o_LIFT;
  occa::memory o_DrT, o_DsT, o_DtT, o_LIFTT;

  occa::memory o_D; // tensor product differentiation matrix (for Hexes)
  
  occa::memory o_vgeo, o_sgeo;
  occa::memory o_vmapM, o_vmapP, o_mapP;
  
  occa::memory o_EToE, o_EToF;
  occa::memory o_EToB, o_x, o_y, o_z;

  // cubature (for wadg)
  occa::memory o_intLIFTT, o_intInterpT, o_intx, o_inty;
  occa::memory *o_cubDrWT, *o_cubDsWT, *o_cubDtWT;
  occa::memory *o_cubInterpT, *o_cubProjectT;
  occa::memory o_invMc; // for comparison: inverses of weighted mass matrices
  occa::memory o_c2;

  //MRAB element lists
  occa::memory *o_MRABelementIds, *o_MRABhaloIds;
  occa::memory *o_MRABpmlElementIds, *o_MRABpmlIds;
  occa::memory *o_MRABpmlHaloElementIds, *o_MRABpmlHaloIds;

  occa::memory *o_MRABNelP, *o_MRABNhaloEleP;
  occa::memory **o_MRABelIdsP, **o_MRABhaloIdsP;

  occa::memory *o_MRABpmlNelP, *o_MRABpmlNhaloEleP;
  occa::memory **o_MRABpmlElIdsP, **o_MRABpmlHaloEleIdsP;
  occa::memory **o_MRABpmlIdsP, **o_MRABpmlHaloIdsP;

  // DG halo exchange info
  occa::memory o_haloElementList;
  occa::memory o_haloBuffer;

  occa::memory o_internalElementIds;
  occa::memory o_notInternalElementIds;
  
  // Bernstein-Bezier occa arrays
  occa::memory *o_D0ids, *o_D1ids, *o_D2ids, *o_D3ids, *o_Dvals; // Bernstein deriv matrix indices
  occa::memory *o_invVB1DT, *o_invVB2DT;
  occa::memory *o_VBq, *o_PBq; // cubature interpolation/projection matrices
  occa::memory *o_L0vals, *o_L0ids, *o_ELids, *o_ELvals; 
  occa::memory *o_BBLower, *o_BBRaiseids, *o_BBRaiseVals; 


  // pml vars
  occa::memory o_sigmax, o_sigmay;

  int pmlNelements;
  int nonPmlNelements;
  occa::memory o_pmlElementIds;
  occa::memory o_nonPmlElementIds;

  occa::memory o_pmlqx, o_rhspmlqx, o_respmlqx;
  occa::memory o_pmlqy, o_rhspmlqy, o_respmlqy;
  occa::memory o_pmlNT, o_rhspmlNT, o_respmlNT;
  
  occa::memory o_pmlElementList;
  occa::memory o_pmlSigmaX, o_pmlSigmaY, o_pmlSigmaZ;
  
  occa::memory o_pmlrhsq, o_pmlresq, o_pmlq;

  

  
  // CG gather-scatter info
  void *gsh; // gslib struct pointer

  int *gatherLocalIds; // local index of rank/gather id sorted list of nodes
  int *gatherBaseIds;  // gather index of ""
  int *gatherBaseRanks; // base rank
  int *gatherMaxRanks;  // maximum rank connected to each sorted node
  int *gatherHaloFlags;  // maximum rank connected to each sorted node

  int NuniqueBases; // number of unique bases on this rank
  occa::memory o_gatherNodeOffsets; // list of offsets into gatherLocalNodes for start of base
  occa::memory o_gatherLocalNodes; // indices of local nodes collected by base node
  occa::memory o_gatherTmp; // temporary array to store base values gathered locally

  int NnodeHalo; // number of halo bases on this rank
  occa::memory o_nodeHaloIds;  // indices of halo base nodes after initial local gather
  occa::memory o_subGatherTmp; // temporary DEVICE array to store halo base values prior to DEVICE>HOST copy
  dfloat        *subGatherTmp; // temporary HALO array

  occa::memory o_ggeo; // second order geometric factors
  occa::memory o_projectL2; // local weights for projection.

  occa::kernel *volumeKernel;
  occa::kernel *surfaceKernel;
  occa::kernel *updateKernel;
  occa::kernel *traceUpdateKernel;
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
  occa::kernel *pmlVolumeKernel;
  occa::kernel *pmlSurfaceKernel;
  occa::kernel *pmlUpdateKernel;
  occa::kernel *pmlTraceUpdateKernel;
  
}mesh_t;

// serial sort
void mysort(int *data, int N, const char *order);

// sort entries in an array in parallel
void parallelSort(int N, void *vv, size_t sz,
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

/* renumber global nodes to remove gaps */
void meshParallelConsecutiveGlobalNumbering(int Nnum, int *globalNumbering);

void meshHaloSetupP(mesh_t *mesh);

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

extern "C"
{
  void *gsParallelGatherScatterSetup(int Ngather, int *gatherIds);
  void gsParallelGatherScatter(void *gsh, void *v, const char *type, const char *op);
  void gsParallelGatherScatterDestroy(void *gsh);

  void * xxtSetup(uint num_local_rows, 
                  void* row_ids,
                  uint nnz, 
                  void*   A_i,
                  void*   A_j,
                  void* A_vals,
                  int null_space,
                  char* inttype,
                  char* floattype);

  int xxtSolve(void* x,
               void* A,
               void* rhs);

  int xxtFree(void* A) ;
}



#endif

