#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mesh2D.h"


typedef struct {
  iint id;
  iint level;
  dfloat weight;

  // 4 for maximum number of vertices per element in 2D
  iint v[4];
  dfloat EX[4], EY[4];

  iint cRank;
  iint cId;
} cElement_t;

typedef struct {
  iint Nelements;
  iint offSet;
} cluster_t;

void meshBuildMRABClusters2D(mesh_t *mesh, iint lev, dfloat *weights, iint *levels,
            iint *Nclusters, cluster_t **clusters, iint *Nelements, cElement_t **newElements);

// geometric partition of clusters of elements in 2D mesh using Morton ordering + parallelSort
dfloat meshClusteredGeometricPartition2D(mesh2D *mesh, iint Nclusters, cluster_t *clusters, 
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

void meshMRABWeightedPartitionTri2D(mesh2D *mesh, dfloat *weights,
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
  meshBuildMRABClusters2D(mesh, -1, weights, levels, &Nclusters, &clusters, &Nelements, &elements);
  meshClusteredGeometricPartition2D(mesh, Nclusters, clusters, &Nelements, &elements);

  //initialize the accepted partition
  iint acceptedNelements = Nelements;
  acceptedPartition = elements;

  for (iint lev = 0; lev<1; lev++) {
    if (rank==0) printf("Clustering level %d...", lev);
    meshBuildMRABClusters2D(mesh, lev, weights, levels, &Nclusters, &clusters, &Nelements, &elements);
    if (rank==0) printf("done.\n");
    dfloat partQuality = meshClusteredGeometricPartition2D(mesh, Nclusters, clusters, &Nelements, &elements);

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
  free(mesh->EToV);
  free(mesh->EX);
  free(mesh->EY);
  
  mesh->Nelements = acceptedNelements;

  mesh->EToV = (iint*) calloc(mesh->Nelements*mesh->Nverts, sizeof(iint));
  mesh->EX = (dfloat*) calloc(mesh->Nelements*mesh->Nverts, sizeof(dfloat));
  mesh->EY = (dfloat*) calloc(mesh->Nelements*mesh->Nverts, sizeof(dfloat));

  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Nverts;++n){
      mesh->EToV[e*mesh->Nverts + n] = acceptedPartition[e].v[n];
      mesh->EX  [e*mesh->Nverts + n] = acceptedPartition[e].EX[n];
      mesh->EY  [e*mesh->Nverts + n] = acceptedPartition[e].EY[n];
    }
  }

  // connect elements using parallel sort
  meshParallelConnect(mesh);

  // print out connectivity statistics
  meshPartitionStatistics(mesh);
  
  // connect elements to boundary faces
  meshConnectBoundary(mesh);

  // compute physical (x,y) locations of the element nodes
  meshPhysicalNodesTri2D(mesh);

  // compute geometric factors
  meshGeometricFactorsTri2D(mesh);

  // set up halo exchange info for MPI (do before connect face nodes)
  meshHaloSetup(mesh);

  // connect face nodes (find trace indices)
  meshConnectFaceNodes2D(mesh);

  // compute surface geofacs
  meshSurfaceGeometricFactorsTri2D(mesh);

  // global nodes
  meshParallelConnectNodes(mesh);

  free(mesh->MRABlevel);
  mesh->MRABlevel = (iint *) calloc(mesh->Nelements+mesh->totalHaloPairs,sizeof(iint));
  for(iint e=0;e<mesh->Nelements;++e) 
    mesh->MRABlevel[e] = acceptedPartition[e].level;
  
  free(acceptedPartition);
}