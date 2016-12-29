#ifndef MESH2D_H 
#define MESH2D_H 1

#include <math.h>
#include <stdlib.h>
#include <occa.hpp>

#if 1
#define iint int
#define dfloat float
#define MPI_IINT MPI_INT
#define MPI_DFLOAT MPI_FLOAT
#define iintFormat "%d"
#define dfloatFormat "%f"
#else
#define iint int
#define dfloat double
#define MPI_IINT MPI_INT
#define MPI_DFLOAT MPI_DOUBLE
#define iintFormat "%d"
#define dfloatFormat "%lf"
#endif

typedef struct {

  iint Nverts, Nfaces;

  iint Nnodes;
  dfloat *EX; // coordinates of vertices for each element
  dfloat *EY;
  
  iint Nelements;
  iint *EToV; // element-to-vertex connectivity
  iint *EToE; // element-to-element connectivity
  iint *EToF; // element-to-(local)face connectivity
  iint *EToP; // element-to-partition/process connectivity
  iint *EToB; // element-to-boundary condition type

  // boundary faces
  iint NboundaryFaces; // number of boundary faces
  iint *boundaryInfo; // list of boundary faces (type, vertex-1, vertex-2)

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

  // volume node info 
  iint N, Np;
  dfloat *r, *s;   // coordinates of local nodes
  dfloat *Dr, *Ds; // collocation differentiation matrices
  dfloat *x, *y;   // coordinates of physical nodes

  // quad specific quantity
  iint Nq;
  
  dfloat *D; // 1D differentiation matrix (for tensor-product)
  dfloat *gllz; // 1D GLL quadrature nodes
  dfloat *gllw; // 1D GLL quadrature weights

  // face node info
  iint Nfp;        // number of nodes per face
  iint *faceNodes; // list of element reference interpolation nodes on element faces
  iint *vmapM;     // list of volume nodes that are face nodes
  iint *vmapP;     // list of volume nodes that are paired with face nodes

  dfloat *LIFT; // lift matrix

  iint   Nsgeo;
  dfloat *sgeo;

  // field info for PDE solver
  iint Nfields;
  dfloat *q;    // solution data array
  dfloat *rhsq; // right hand side data array
  dfloat *resq; // residual data array (for LSERK time-stepping)
  
  dfloat Lambda2; // square of penalty paramater used in constructing q^*

  // time stepping info
  dfloat dt; // time step
  dfloat finalTime; // final time to run acoustics to
  iint   NtimeSteps;// number of time steps 
  iint   errorStep; // number of steps between error calculations
  iint   Nrk;
  dfloat rka[5], rkb[5], rkc[5];

  // ploting info for generating field vtu
  iint    plotNverts;    // number of vertices for each plot element
  iint    plotNp;        // number of plot nodes per element
  iint    plotNelements; // number of "plot elements" per element
  iint   *plotEToV;      // triangulation of plot nodes
  dfloat *plotR, *plotS; // coordinates of plot nodes in reference element
  dfloat *plotInterp;    // warp & blend to plot node interpolation matrix

  // volume cubature node info
  iint    cubNp; // number of cubature nodes
  dfloat *cubInterp; // interpolate from W&B to cubature nodes
  dfloat *cubDrW;    // 'r' weak differentiation matrix
  dfloat *cubDsW;    // 's' weak differentiation matrix
  dfloat *cubProject; // projection matrix from cubature nodes to W&B nodes
  
  // surface integration node info
  iint    intNfp;    // number of integration nodes on each face
  dfloat *intInterp; // interp from surface node to integration nodes
  dfloat *intLIFT;   // lift from surface integration nodes to W&B volume nodes
  dfloat *intx, *inty; // coordinates of suface integration nodes
  
  // global numbering info
  iint   *baseIds;   // local index of base nodes for each interp node
  iint   *baseRanks; // rank of base node for each interp node
  
  // elliptic solver info
  iint  NgatherNodes;
  iint  *globalNumbering;
  iint  *gatherCounts;
  iint  *gatherOffsets;
  iint  *gatherIds;

  // elliptic gather info
  iint gatherNsend;
  iint gatherNrecv;
  iint gatherNbaseRecv;

  iint *gatherRecvIds;
  iint *gatherSendIds;  
  iint *gatherRecvCounts;
  iint *gatherSendCounts;
  iint *gatherRecvDispls;
  iint *gatherSendDispls;

  iint *gatherRecvStarts;
  iint *gatherRecvSourceIds;
  iint *gatherRecvBaseIds;

  iint *gatherLocalStarts;
  iint *gatherLocalSourceIds;
  iint *gatherLocalBaseIds;

  // elliptic scatter info
  iint scatterNsend;
  iint scatterNrecv;
  iint scatterNbaseSend;

  iint *scatterRecvIds;
  iint *scatterSendIds;  
  iint *scatterRecvCounts;
  iint *scatterSendCounts;
  iint *scatterRecvDispls;
  iint *scatterSendDispls;

  iint *scatterSendStarts;
  iint *scatterSendSourceIds;
  iint *scatterSendBaseIds;

  iint *scatterLocalStarts;
  iint *scatterLocalSourceIds;
  iint *scatterLocalBaseIds;

  // pml stuff
  iint    pmlNfields;
  iint    pmlNelements;
  iint   *pmlElementList;
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

  // Boltzmann specific stuff
  dfloat RT, sqrtRT, tauInv;
  
  // occa stuff
  occa::device device;
  occa::memory o_q, o_rhsq, o_resq;

  occa::memory o_Dr, o_Ds, o_LIFT;
  occa::memory o_DrT, o_DsT, o_D, o_LIFTT;

  occa::memory o_intLIFTT, o_intInterpT, o_intx, o_inty;
  occa::memory o_cubDrWT, o_cubDsWT, o_cubInterpT, o_cubProjectT;
  
  occa::memory o_vgeo, o_sgeo;
  occa::memory o_vmapM, o_vmapP;

  occa::memory o_EToB, o_x, o_y;
  
  occa::memory o_haloElementList;
  occa::memory o_haloBuffer;

  occa::memory o_internalElementIds;
  occa::memory o_notInternalElementIds;

  // pml vars
  occa::memory o_sigmax, o_sigmay;
  occa::memory o_pmlqx, o_rhspmlqx, o_respmlqx;
  occa::memory o_pmlqy, o_rhspmlqy, o_respmlqy;
  occa::memory o_pmlNT, o_rhspmlNT, o_respmlNT;
  
  occa::memory o_pmlElementList;
  occa::memory o_pmlSigmaX, o_pmlSigmaY;
  
  occa::memory o_pmlrhsq, o_pmlresq, o_pmlq;

  occa::kernel haloExtractKernel;
  
  occa::kernel volumeKernel;
  occa::kernel surfaceKernel;
  occa::kernel partialSurfaceKernel;
  occa::kernel updateKernel;
  occa::kernel relaxationKernel;

  occa::kernel pmlKernel;
  occa::kernel pmlUpdateKernel;
  
}mesh2D;

mesh2D* meshReader2D(char *fileName);

// mesh readers
mesh2D* meshParallelReader2D(char *fileName);
mesh2D* meshParallelReaderQuad2D(char *fileName);

// build connectivity in serial
void meshConnect2D(mesh2D *mesh);

// build element-boundary connectivity
void meshConnectBoundary2D(mesh2D *mesh);

// build connectivity in parallel
void meshParallelConnect2D(mesh2D *mesh);

// create global number of nodes
void meshNumberNodes2D(mesh2D *mesh);

// repartition elements in parallel
void meshGeometricPartition2D(mesh2D *mesh);

// print out mesh 
void meshPrint2D(mesh2D *mesh);

// print out mesh in parallel from the root process
void meshParallelPrint2D(mesh2D *mesh);

// print out mesh partition in parallel
void meshVTU2D(mesh2D *mesh, char *fileName);

// print out solution at plot nodes 
void meshPlotVTU2D(mesh2D *mesh, char *fileNameBase, iint fld);

// sort entries in an array in parallel
void parallelSort(iint N, void *vv, size_t sz,
		  int (*compare)(const void *, const void *),
		  void (*match)(void *, void *)
		  );

// compute geometric factors for local to physical map 
void meshGeometricFactors2D(mesh2D *mesh);
void meshGeometricFactorsQuad2D(mesh2D *mesh);

void meshSurfaceGeometricFactors2D(mesh2D *mesh);
void meshSurfaceGeometricFactorsQuad2D(mesh2D *mesh);

void meshPhysicalNodes2D(mesh2D *mesh);
void meshPhysicalNodesQuad2D(mesh2D *mesh);

void meshLoadReferenceNodes2D(mesh2D *mesh, int N);
void meshLoadReferenceNodesQuad2D(mesh2D *mesh, int N);

void meshGradient2D(mesh2D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy);
void meshGradientQuad2D(mesh2D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy);

// print out parallel partition i
void meshPartitionStatistics2D(mesh2D *mesh);

// functions that call OCCA kernels
void occaTest(mesh2D *mesh);

// 
void occaOptimizeGradient2D(mesh2D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy);
void occaOptimizeGradientQuad2D(mesh2D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy);

// serial face-node to face-node connection
void meshConnectFaceNodes2D(mesh2D *mesh);

// halo connectivity information
void meshHaloSetup2D(mesh2D *mesh);

// perform complete halo exchange
void meshHaloExchange2D(mesh2D *mesh,
			size_t Nbytes, // number of bytes per element
			void *sourceBuffer, 
			void *sendBuffer, 
			void *recvBuffer);

// start halo exchange
void meshHaloExchangeStart2D(mesh2D *mesh,
			     size_t Nbytes,       // message size per element
			     void *sendBuffer,    // outgoing halo
			     void *recvBuffer);   // incoming halo

// finish halo exchange
void meshHaloExchangeFinish2D(mesh2D *mesh);

// extract halo data from sourceBuffer and save to sendBuffer
void meshHaloExtract2D(mesh2D *mesh, size_t Nbytes, void *sourceBuffer, void *sendBuffer);

// build list of nodes on each face of the reference element
void meshBuildFaceNodes2D(mesh2D *mesh);
void meshBuildFaceNodesQuad2D(mesh2D *mesh);

mesh2D *meshSetup2D(char *filename, int N);
mesh2D *meshSetupQuad2D(char *filename, int N);

void meshAcousticsRun2D(mesh2D *mesh);
void meshAcousticsRunQuad2D(mesh2D *mesh);
void meshAcousticsOccaRun2D(mesh2D *mesh);
void meshAcousticsOccaAsyncRun2D(mesh2D *mesh);
void meshAcousticsSplitSurfaceOccaAsyncRun2D(mesh2D *mesh);

void meshAcousticsSetup2D(mesh2D *mesh);
void meshAcousticsSetupQuad2D(mesh2D *mesh);
void meshAcousticsVolume2D(mesh2D *mesh);
void meshAcousticsVolumeQuad2D(mesh2D *mesh);
void meshAcousticsSurface2D(mesh2D *mesh, dfloat time);
void meshAcousticsSurfaceQuad2D(mesh2D *mesh);
void meshAcousticsUpdate2D(mesh2D *mesh, dfloat rka, dfloat rkb);
void meshAcousticsError2D(mesh2D *mesh, dfloat time);

// set up perfectly matched layer
void meshAcousticsPmlSetup2D(mesh2D *mesh,
			     dfloat xmin, dfloat xmax, // bounding box for non-pml sub-domain
			     dfloat ymin, dfloat ymax,
			     dfloat xsigma, dfloat ysigma);

// add PML right hand stuff
void meshAcousticsPml2D(mesh2D *mesh);
void meshAcousticsPmlUpdate2D(mesh2D *mesh, dfloat rka, dfloat rkb);

// compute solution to cavity problem
void acousticsCavitySolution2D(dfloat x, dfloat y, dfloat time, 
			       dfloat *u, dfloat *v, dfloat *p);

// initial Gaussian pulse
void acousticsGaussianPulse2D(dfloat x, dfloat y, dfloat t,
			      dfloat *u, dfloat *v, dfloat *p);

void meshAcousticsComputeVorticity2D(mesh2D *mesh, dfloat *q, iint outfld, iint Nfields);

// Boltzmann model
void meshBoltzmannSetup2D(mesh2D *mesh);
void meshBoltzmannVolume2D(mesh2D *mesh);
void meshBoltzmannSurface2D(mesh2D *mesh, dfloat t);
void meshBoltzmannUpdate2D(mesh2D *mesh, dfloat rka, dfloat rkb);
void meshBoltzmannRun2D(mesh2D *mesh);
void meshBoltzmannOccaRun2D(mesh2D *mesh);
void meshBoltzmannError2D(mesh2D *mesh, dfloat time);

void boltzmannGaussianPulse2D(dfloat x, dfloat y, dfloat t,
			      dfloat *q1, dfloat *q2, dfloat *q3,
			      dfloat *q4, dfloat *q5, dfloat *q6);

void meshBoltzmannComputeVorticity2D(mesh2D *mesh, dfloat *q, iint outfld, iint Nfields);

// continuous Galerkin functions
void meshGatherScatterSetupQuad2D(mesh2D *mesh);
void meshGatherScatterQuad2D(mesh2D *mesh, dfloat *AqL);
void meshEllipticSetupQuad2D(mesh2D *mesh);
void meshEllipticSolveQuad2D(mesh2D *mesh);
void meshEllipticLocalOpQuad2D(mesh2D *mesh, dfloat *qL, dfloat lambda, dfloat *AqL);

#define mymax(a,b) ((a>b)?a:b)
#define mymin(a,b) ((a<b)?a:b)
#define norm(a,b) ( sqrt((a)*(a)+(b)*(b)) )

/* offsets for geometric factors */
#define RXID 0  
#define RYID 1  
#define SXID 2  
#define SYID 3  
#define  JID 4
#define JWID 5

/* offsets for nx, ny, sJ, 1/J */
#define NXID 0  
#define NYID 1  
#define SJID 2  
#define IJID 3  


#endif

