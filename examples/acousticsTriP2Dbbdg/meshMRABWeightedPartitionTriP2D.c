#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mesh2D.h"


typedef struct {
  iint id;
  iint level;
  dfloat weight;
  iint N;

  // 4 for maximum number of vertices per element in 2D
  iint v[4];
  dfloat EX[4], EY[4];

  iint cRank;
  iint cId;
  int type;
} cElement_t;

typedef struct {
  iint Nelements;
  iint offSet;
} cluster_t;

void meshBuildMRABClustersP2D(mesh_t *mesh, iint lev, dfloat *weights, iint *levels,
            iint *Nclusters, cluster_t **clusters, iint *Nelements, cElement_t **newElements);

// geometric partition of clusters of elements in 2D mesh using Morton ordering + parallelSort
dfloat meshClusteredGeometricPartitionP2D(mesh2D *mesh, iint Nclusters, cluster_t *clusters, 
                              iint *Nelements, cElement_t **elements);

/* ---------------------------------------------------------

This function is a bit spaghetti, but the general idea is 
we cluster low-MRAB-level elements together along with a 
halo and partition the mesh of clusters. This reduces the MPI 
costs of communicating on the low levels.

The algorithm performs the following steps
  - cluster elements of level lev or lower
  - put clusters together a single 'owning' process
  - sort the list of clusters using a space-filling curve
  - partition the SFC between the processors, exchange the
    elements along the processor boundaries to improve the
    partitioning.
  - If the resulting partition is acceptable, save it.
  - If not, return to the last acceptable partition, and rerun 
    the mesh setup. 

------------------------------------------------------------ */

void meshMRABWeightedPartitionTriP2D(mesh2D *mesh, dfloat *weights,
                                      iint numLevels, iint *levels) {

  const dfloat TOL = 0.8; //tolerance on what partitions are ruled 'acceptable'
                          // min_{ranks} totalWeight > TOL*max_{ranks} totalWeight => accepted

  iint rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  iint Nelements, Nclusters;

  cElement_t *elements, *acceptedPartition;
  cluster_t *clusters;

  if (!levels) numLevels = 1;

  //perform the first weigthed partitioning with no clustering
  meshBuildMRABClustersP2D(mesh, -1, weights, levels, &Nclusters, &clusters, &Nelements, &elements);
  meshClusteredGeometricPartitionP2D(mesh, Nclusters, clusters, &Nelements, &elements);

  //initialize the accepted partition
  iint acceptedNelements = Nelements;
  acceptedPartition = elements;

  for (iint lev = 0; lev<mesh->MRABNlevels; lev++) {
    if (rank==0) printf("Clustering level %d...", lev);
    meshBuildMRABClustersP2D(mesh, lev, weights, levels, &Nclusters, &clusters, &Nelements, &elements);
    if (rank==0) printf("done.\n");
    dfloat partQuality = meshClusteredGeometricPartitionP2D(mesh, Nclusters, clusters, &Nelements, &elements);

    if (partQuality > TOL) {
      if (rank ==0) printf("Accepting level %d clustered partition...(quality = %g)\n", lev, partQuality);
      free(acceptedPartition); //discard the old partition
      acceptedNelements = Nelements;
      acceptedPartition = elements; //good partition
    } else {
      if (rank ==0) printf("Regecting level %d clustered partition...(quality = %g)\n", lev, partQuality);
      free(elements); //discard this partition
      break;  
    }
  }

  //save this partition, and perform the mesh setup again.   
  mesh->Nelements = acceptedNelements;

  mesh->EToV = (iint*) realloc(mesh->EToV, mesh->Nelements*mesh->Nverts*sizeof(iint));
  mesh->EX = (dfloat*) realloc(mesh->EX, mesh->Nelements*mesh->Nverts*sizeof(dfloat));
  mesh->EY = (dfloat*) realloc(mesh->EY, mesh->Nelements*mesh->Nverts*sizeof(dfloat));
  mesh->N  =   (iint*) realloc(mesh->N,  mesh->Nelements*sizeof(iint));
  mesh->elementInfo = (int *) realloc(mesh->elementInfo,mesh->Nelements*sizeof(int));
  mesh->MRABlevel = (iint *) realloc(mesh->MRABlevel,mesh->Nelements*sizeof(iint));

  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Nverts;++n){
      mesh->EToV[e*mesh->Nverts + n] = acceptedPartition[e].v[n];
      mesh->EX  [e*mesh->Nverts + n] = acceptedPartition[e].EX[n];
      mesh->EY  [e*mesh->Nverts + n] = acceptedPartition[e].EY[n];
    }
    mesh->N[e] = acceptedPartition[e].N;
    mesh->elementInfo[e] = acceptedPartition[e].type;
    mesh->MRABlevel[e] = acceptedPartition[e].level;
  }

  // connect elements using parallel sort
  meshParallelConnect(mesh);

  // print out connectivity statistics
  meshPartitionStatistics(mesh);
  
  // connect elements to boundary faces
  meshConnectBoundary(mesh);

  // compute geometric factors
  meshGeometricFactorsTri2D(mesh);

  // compute surface geofacs
  meshSurfaceGeometricFactorsTriP2D(mesh);  

  // set up halo exchange info for MPI (do before connect face nodes)
  meshHaloSetupP(mesh);
  
  // compute physical (x,y) locations of the element nodes
  meshPhysicalNodesTriP2D(mesh);

  // connect face nodes (find trace indices)
  meshConnectFaceNodesP2D(mesh);

  if (mesh->totalHaloPairs) {
    mesh->N  =   (iint*) realloc(mesh->N, (mesh->Nelements+mesh->totalHaloPairs)*sizeof(iint));
    mesh->MRABlevel = (iint *) realloc(mesh->MRABlevel,(mesh->Nelements+mesh->totalHaloPairs)*sizeof(iint));
    iint *MRABsendBuffer = (iint *) calloc(mesh->totalHaloPairs,sizeof(iint));
    meshHaloExchange(mesh, sizeof(iint), mesh->MRABlevel, MRABsendBuffer, mesh->MRABlevel+mesh->Nelements);
    meshHaloExchange(mesh, sizeof(iint), mesh->N, MRABsendBuffer, mesh->N+mesh->Nelements);
    free(MRABsendBuffer);
  }

  
  free(acceptedPartition);
}