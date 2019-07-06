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

#include "mesh.hpp"

typedef struct {
  int id;
  int level;
  dfloat weight;

  // 8 for maximum number of vertices per element in 3D
  int v[8];
  dfloat EX[8], EY[8], EZ[8];

  int cRank;
  int cId;
  int type;
} cElement_t;

typedef struct {
  int Nelements;
  int offSet;
} cluster_t;


void meshBuildMRABClusters2D(mesh_t *mesh, int lev, dfloat *weights, int *levels,
            int *Nclusters, cluster_t **clusters, int *newNelements, cElement_t **elements);

void meshBuildMRABClusters3D(mesh_t *mesh, int lev, dfloat *weights, int *levels,
            int *Nclusters, cluster_t **clusters, int *newNelements, cElement_t **elements);

dfloat meshClusteredGeometricPartition2D(mesh_t *mesh, int Nclusters, cluster_t *clusters,
           int *Nelements, cElement_t **elements);

dfloat meshClusteredGeometricPartition3D(mesh_t *mesh, int Nclusters, cluster_t *clusters,
           int *Nelements, cElement_t **elements);

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
void mesh_t::MRABWeightedPartition(dfloat *weights,
                                   int numLevels, int *levels) {

  const dfloat TOL = 0.8; //tolerance on what partitions are ruled 'acceptable'
                          // min_{ranks} totalWeight > TOL*max_{ranks} totalWeight => accepted

  int newNelements, Nclusters;

  cElement_t *elements, *acceptedPartition;
  cluster_t *clusters;

  if (!levels) numLevels = 1;

  //perform the first weigthed partitioning with no clustering
  if (dim==2) {
    meshBuildMRABClusters2D(this, -1, weights, levels, &Nclusters, &clusters, &newNelements, &elements);
    meshClusteredGeometricPartition2D(this, Nclusters, clusters, &newNelements, &elements);
  } else {
    meshBuildMRABClusters3D(this, -1, weights, levels, &Nclusters, &clusters, &newNelements, &elements);
    meshClusteredGeometricPartition3D(this, Nclusters, clusters, &newNelements, &elements);
  }

  //initialize the accepted partition
  int acceptedNelements = newNelements;
  acceptedPartition = elements;

  for (int lev = 0; lev<MRABNlevels; lev++) {
    if (rank==0) printf("Clustering level %d...", lev);

    if (dim==2)
      meshBuildMRABClusters2D(this,lev, weights, levels, &Nclusters, &clusters, &newNelements, &elements);
    else
      meshBuildMRABClusters3D(this,lev, weights, levels, &Nclusters, &clusters, &newNelements, &elements);

    if (rank==0) printf("done.\n");
    dfloat partQuality;
    if (dim==2)
      partQuality = meshClusteredGeometricPartition2D(this, Nclusters, clusters, &newNelements, &elements);
    else
      partQuality = meshClusteredGeometricPartition3D(this, Nclusters, clusters, &newNelements, &elements);

    if (partQuality > TOL) {
      if (rank ==0) printf("Accepting level %d clustered partition...(quality = %g)\n", lev, partQuality);
      free(acceptedPartition); //discard the old partition
      acceptedNelements = newNelements;
      acceptedPartition = elements; //good partition
    } else {
      if (rank ==0) printf("Regecting level %d clustered partition...(quality = %g)\n", lev, partQuality);
      free(elements); //discard this partition
      break;
    }
  }

  //save this partition, and perform the mesh setup again.
  Nelements = acceptedNelements;

  EToV = (hlong*) realloc(EToV,Nelements*Nverts*sizeof(hlong));
  EX = (dfloat*) realloc(EX,Nelements*Nverts*sizeof(dfloat));
  EY = (dfloat*) realloc(EY,Nelements*Nverts*sizeof(dfloat));
  if (dim==3)
    EZ = (dfloat*) realloc(EZ,Nelements*Nverts*sizeof(dfloat));

  elementInfo = (hlong *) realloc(elementInfo,Nelements*sizeof(hlong));
  MRABlevel = (int *) realloc(MRABlevel,Nelements*sizeof(int));

  for(dlong e=0;e<Nelements;++e){
    for(int n=0;n<Nverts;++n){
      EToV[e*Nverts + n] = acceptedPartition[e].v[n];
      EX  [e*Nverts + n] = acceptedPartition[e].EX[n];
      EY  [e*Nverts + n] = acceptedPartition[e].EY[n];
      if (dim==3)
        EZ  [e*Nverts + n] = acceptedPartition[e].EZ[n];
    }
    elementInfo[e] = acceptedPartition[e].type;
    MRABlevel[e] = acceptedPartition[e].level;
  }

  // connect elements using parallel sort
  this->ParallelConnect();

  // print out connectivity statistics
  this->PrintPartitionStatistics();

  // connect elements to boundary faces
  this->ConnectBoundary();


  if(dim==3 && NfaceVertices==2) // Tri/Quad 3D
    this->LoadReferenceNodes(N);

  this->PhysicalNodes();
  this->GeometricFactors();

  // set up halo exchange info for MPI (do before connect face nodes)
  this->HaloSetup();

  // connect face nodes (find trace indices)
  this->ConnectFaceNodes();

  this->SurfaceGeometricFactors();

  // global nodes
  this->ParallelConnectNodes();

  if (totalHaloPairs) {
    MRABlevel = (int *) realloc(MRABlevel,(Nelements+totalHaloPairs)*sizeof(int));
    halo->Exchange(MRABlevel, 1, ogs_int);
  }

  free(acceptedPartition);
}

