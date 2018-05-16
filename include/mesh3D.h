#ifndef MESH3D_H 
#define MESH3D_H 1

// generic mesh structure 
#include "mesh.h"

#define mesh3D mesh_t

// mesh readers
mesh3D* meshParallelReaderTri3D(char *fileName);
mesh3D* meshParallelReaderQuad3D(char *fileName);
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

// compute geometric factors for local to physical map
void meshGeometricFactorsTri3D(mesh3D *mesh);
void meshGeometricFactorsQuad3D(mesh3D *mesh);
void meshGeometricFactorsTet3D(mesh3D *mesh);
void meshGeometricFactorsHex3D(mesh3D *mesh);

void meshSurfaceGeometricFactorsTri3D(mesh3D *mesh);
void meshSurfaceGeometricFactorsQuad3D(mesh3D *mesh);
void meshSurfaceGeometricFactorsTet3D(mesh3D *mesh);
void meshSurfaceGeometricFactorsHex3D(mesh3D *mesh);

void meshPhysicalNodesTri3D(mesh3D *mesh);
void meshPhysicalNodesQuad3D(mesh3D *mesh);
void meshSphericalNodesQuad3D(mesh3D *mesh);
void meshEquiSphericalExtensionQuad3D(mesh3D *mesh);
void meshExtendGridQuad3D(mesh3D *mesh);
void meshPreserveGridQuad3D(mesh3D *mesh);
void meshEquiSphericalNodesQuad3D(mesh3D *mesh);
void meshPhysicalNodesTet3D(mesh3D *mesh);
void meshPhysicalNodesHex3D(mesh3D *mesh);

void meshLoadReferenceNodesTet3D(mesh3D *mesh, int N);
void meshLoadReferenceNodesHex3D(mesh3D *mesh, int N);

void meshGradientTet3D(mesh3D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy, dfloat *dqdz);
void meshGradientHex3D(mesh3D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy, dfloat *dqdz);

// print out parallel partition i
void meshPartitionStatistics3D(mesh3D *mesh);

// default occa set up
void meshOccaSetup3D(mesh3D *mesh, char *deviceConfig, occa::kernelInfo &kernelInfo);

// functions that call OCCA kernels
void occaTest3D(mesh3D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy, dfloat *dqdz);

// 
void occaOptimizeGradientTet3D(mesh3D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy, dfloat *dqdz);
void occaOptimizeGradientHex3D(mesh3D *mesh, dfloat *q, dfloat *dqdx, dfloat *dqdy, dfloat *dqdz);

// serial face-node to face-node connection
void meshConnectFaceNodes3D(mesh3D *mesh);

//
mesh3D *meshSetupTri3D(char *filename, int N, dfloat sphereRadius);
mesh3D *meshSetupQuad3D(int mesh_size, int N, dfloat sphereRadius, char *mode);
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

void meshParallelGatherScatter(mesh3D *mesh,
			       ogs_t *ogs, 
			       occa::memory &o_v,
			       occa::memory &o_gsv,
			       const char *type,
			       const char *op);

ogs_t *meshParallelGatherScatterSetup(mesh3D *mesh,    // provides DEVICE
				      iint Nlocal,     // number of local nodes
				      iint Nbytes,     // number of bytes per node
				      iint *localIds,  // local index of nodes
				      iint *baseIds,   // gather index of their base nodes
				      iint *haloFlags); // 1 for halo node, 0 for not

void meshMRABSetup3D(mesh3D *mesh, dfloat *EToDT, int maxLevels);
void meshMRABSetupQuad3D(mesh3D *mesh, dfloat *EToDT, int maxLevels);

//MRAB weighted mesh partitioning
void meshMRABWeightedPartitionTet3D(mesh3D *mesh, dfloat *weights,
                                      iint numLevels, iint *levels);

void meshMRABWeightedPartitionQuad3D(mesh3D *mesh, dfloat *weights,
				     iint numLevels, iint *levels);

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
#define IHID 5
#define WSJID 6
//
//offsets for boltzmann PML variables
#define QXID1 0  
#define QXID2 1  
#define QXID3 2
#define QXID4 3  
#define QXID5 4  
#define QXID6 5  
#define QXID8 6 
//
#define QYID1 7  
#define QYID2 8  
#define QYID3 9
#define QYID4 10  
#define QYID5 11  
#define QYID7 12  
#define QYID9 13 
//
#define QZID1 14  
#define QZID2 15  
#define QZID3 16
#define QZID4 17  
#define QZID6 18  
#define QZID7 19  
#define QZID10  20   

#endif

