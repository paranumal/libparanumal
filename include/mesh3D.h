#ifndef MESH3D_H 
#define MESH3D_H 1

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

// OCCA+gslib gather scatter
typedef struct {
  
  iint         Ngather;          //  number of gather nodes

  iint         *gatherOffsets;
  iint         *gatherHaloFlags;
  iint         *gatherBaseRanks;
  iint         *gatherLocalIds;
  iint         *gatherBaseIds;
  
  occa::memory o_gatherOffsets;  //  start of local bases
  occa::memory o_gatherLocalIds; //  base connected nodes
  occa::memory o_gatherTmp;      //  DEVICE gather buffer
  void         *gatherGsh;       // gslib gather 

  iint         Nscatter;
  occa::memory o_scatterOffsets; //  start of local bases
  occa::memory o_scatterLocalIds;//  base connected nodes
  
  iint         Nhalo;            //  number of halo nodes
  occa::memory o_haloLocalIds;   //  list of halo nodes to
  occa::memory o_haloTmp;        //  temporary halo buffer
  void         *haloTmp;         //  temporary HOST halo buffer

}ogs_t;

typedef struct {

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
  
  // face node info
  iint Nfp;        // number of nodes per face
  iint *faceNodes; // list of element reference interpolation nodes on element faces
  iint *vmapM;     // list of volume nodes that are face nodes
  iint *vmapP;     // list of volume nodes that are paired with face nodes
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

  // overlapping additive schwarz direction transform (Hex or quad only)
  
  
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
  occa::memory o_cubInterpT, o_cubProjectT;
  occa::memory o_invMc; // for comparison: inverses of weighted mass matrices
  occa::memory o_c2;

  // DG halo exchange info
  occa::memory o_haloElementList;
  occa::memory o_haloBuffer;

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

  occa::kernel gatherKernel;
  occa::kernel scatterKernel;

  occa::kernel getKernel;
  occa::kernel putKernel;

  occa::kernel AxKernel;
  occa::kernel weightedInnerProduct1Kernel;
  occa::kernel weightedInnerProduct2Kernel;
  occa::kernel scaledAddKernel;
  occa::kernel dotMultiplyKernel;
  occa::kernel dotDivideKernel;
}mesh3D;

// mesh readers
mesh3D* meshParallelReaderTet3D(char *fileName);
mesh3D* meshParallelReaderHex3D(char *fileName);

// build connectivity in serial
void meshConnect3D(mesh3D *mesh);

// build element-boundary connectivity
void meshConnectBoundary3D(mesh3D *mesh);

// build connectivity in parallel
void meshParallelConnect3D(mesh3D *mesh);

// repartition elements in parallel
void meshGeometricPartition3D(mesh3D *mesh);

// print out mesh 
void meshPrint3D(mesh3D *mesh);

// print out mesh in parallel from the root process
void meshParallelPrint3D(mesh3D *mesh);

// print out mesh partition in parallel
void meshVTU3D(mesh3D *mesh, char *fileName);

// print out mesh field
void meshPlotVTU3D(mesh3D *mesh, char *fileNameBase, iint fld);

// serial sort
void mysort(iint *data, iint N, const char *order);

// sort entries in an array in parallel
void parallelSort(iint N, void *vv, size_t sz,
		  int (*compare)(const void *, const void *),
		  void (*match)(void *, void *)
		  );

// compute geometric factors for local to physical map 
void meshGeometricFactorsTet3D(mesh3D *mesh);
void meshGeometricFactorsHex3D(mesh3D *mesh);

void meshSurfaceGeometricFactorsTet3D(mesh3D *mesh);
void meshSurfaceGeometricFactorsHex3D(mesh3D *mesh);

void meshPhysicalNodesTet3D(mesh3D *mesh);
void meshPhysicalNodesHex3D(mesh3D *mesh);

void meshLoadReferenceNodesTet3D(mesh3D *mesh, int N);
void meshLoadReferenceNodesHex3D(mesh3D *mesh, int N);

void meshGradientTet3D(mesh3D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy, dfloat *dqdz);
void meshGradientHex3D(mesh3D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy, dfloat *dqdz);

// print out parallel partition i
void meshPartitionStatistics3D(mesh3D *mesh);

// functions that call OCCA kernels
void occaTest3D(mesh3D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy, dfloat *dqdz);

// 
void occaOptimizeGradientTet3D(mesh3D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy, dfloat *dqdz);
void occaOptimizeGradientHex3D(mesh3D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy, dfloat *dqdz);

// serial face-node to face-node connection
void meshConnectFaceNodes3D(mesh3D *mesh);

//
mesh3D *meshSetupTet3D(char *filename, int N);
mesh3D *meshSetupHex3D(char *filename, int N);

void meshParallelConnectNodesHex3D(mesh3D *mesh);


// halo connectivity information
void meshHaloSetup3D(mesh3D *mesh);

// perform halo exchange
void meshHaloExchange3D(mesh3D *mesh,
			size_t Nbytes,  // number of bytes per element
			void *sourceBuffer, 
			void *sendBuffer, 
			void *recvBuffer);

void meshHaloExchangeStart3D(mesh3D *mesh,
			     size_t Nbytes,       // message size per element
			     void *sendBuffer,    // temporary buffer
			     void *recvBuffer);

void meshHaloExchangeFinish3D(mesh3D *mesh);

// build list of nodes on each face of the reference element
void meshBuildFaceNodes3D(mesh3D *mesh);
void meshBuildFaceNodesHex3D(mesh3D *mesh);



// void meshParallelGatherScatter3D(mesh3D *mesh, occa::memory &o_v, occa::memory &o_gsv, const char *type);

void meshParallelGatherScatter3D(mesh3D *mesh,
				 ogs_t *ogs, 
				 occa::memory &o_v,
				 occa::memory &o_gsv,
				 const char *type);

ogs_t *meshParallelGatherScatterSetup3D(mesh3D *mesh,    // provides DEVICE
					iint Nlocal,     // number of local nodes
					iint Nbytes,     // number of bytes per node
					iint *localIds,  // local index of nodes
					iint *baseIds,   // gather index of their base nodes
					iint *baseRanks, // rank of their base nodes
					iint *haloFlags); // 1 for halo node, 0 for not

void meshAcousticsRun3D(mesh3D *mesh);
void meshAcousticsSetup3D(mesh3D *mesh);
void meshAcousticsVolume3D(mesh3D *mesh);
void meshAcousticsSurface3D(mesh3D *mesh, dfloat time);
void meshAcousticsUpdate3D(mesh3D *mesh, dfloat rka, dfloat rkb);
void meshAcousticsError3D(mesh3D *mesh, dfloat time);

void meshAcousticsOccaRun3D(mesh3D *mesh);
void meshAcousticsOccaRunHex3D(mesh3D *mesh);

void acousticsCavitySolution3D(dfloat x, dfloat y, dfloat z, dfloat time,
			       dfloat *u, dfloat *v, dfloat *w, dfloat *p);

#define mymax(a,b) (((a)>(b))?(a):(b))
#define mymin(a,b) (((a)<(b))?(a):(b))
#define norm(a,b,c) ( sqrt((a)*(a)+(b)*(b)+(c)*(c)) )

/* offsets for geometric factors */
#define RXID 0  
#define RYID 1  
#define RZID 2
#define SXID 3  
#define SYID 4  
#define SZID 5  
#define TXID 6  
#define TYID 7  
#define TZID 8  
#define  JID 9
#define JWID 10

/* offsets for second order geometric factors */
#define G00ID 0  
#define G01ID 1  
#define G02ID 2
#define G11ID 3  
#define G12ID 4  
#define G22ID 5  
#define GWJID 6  

/* offsets for nx, ny, sJ, 1/J */
#define NXID 0  
#define NYID 1  
#define NZID 2 
#define SJID 3  
#define IJID 4  
#endif

