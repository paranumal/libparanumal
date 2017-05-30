
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mesh3D.h"

#define bitRange 10

typedef struct {
  iint id;
  iint level;
  dfloat weight;

  // 8 for maximum number of vertices per element in 3D
  iint v[8];
  dfloat EX[8], EY[8], EZ[8];

  iint cRank;
  iint cId;
} cElement_t;

typedef struct {
  iint Nelements;
  iint offSet;
} cluster_t;

typedef struct {
  iint Nelements;
  iint offSet;
  iint rank;

  iint destId;
  iint destOffset;
  iint destRank;

  dfloat weight;
  unsigned long long int index; //morton index
} parallelCluster_t;

//This is linked form meshGeometricPartition2D.c
unsigned long long int mortonIndex3D(unsigned int ix, unsigned int iy, unsigned int iz);
void bogusMatch(void *a, void *b);

void meshBuildMRABClusters(mesh_t *mesh, iint lev, dfloat *weights, iint *levels,
            iint *Nclusters, cluster_t **clusters, iint *Nelements, cElement_t **newElements);

dfloat improveClusteredPartition(iint *Nclusters, parallelCluster_t **parallelClusters);

// geometric partition of clusters of elements in 3D mesh using Morton ordering + parallelSort
dfloat meshClusteredGeometricPartition3D(mesh3D *mesh, iint Nclusters, cluster_t *clusters, 
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
void meshMRABWeightedPartition3D(mesh3D *mesh, dfloat *weights,
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
  meshBuildMRABClusters(mesh, -1, weights, levels, &Nclusters, &clusters, &Nelements, &elements);
  meshClusteredGeometricPartition3D(mesh, Nclusters, clusters, &Nelements, &elements);

  //initialize the accepted partition
  iint acceptedNelements = Nelements;
  acceptedPartition = elements;

  for (iint lev = 0; lev<1; lev++) {
    if (rank==0) printf("Clustering level %d...", lev);
    meshBuildMRABClusters(mesh, lev, weights, levels, &Nclusters, &clusters, &Nelements, &elements);
    if (rank==0) printf("done.\n");
    dfloat partQuality = meshClusteredGeometricPartition3D(mesh, Nclusters, clusters, &Nelements, &elements);

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

  mesh->EToV = (iint*) realloc(mesh->EToV,mesh->Nelements*mesh->Nverts*sizeof(iint));
  mesh->EX = (dfloat*) realloc(mesh->EX,mesh->Nelements*mesh->Nverts*sizeof(dfloat));
  mesh->EY = (dfloat*) realloc(mesh->EY,mesh->Nelements*mesh->Nverts*sizeof(dfloat));
  mesh->EZ = (dfloat*) realloc(mesh->EZ,mesh->Nelements*mesh->Nverts*sizeof(dfloat));

  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Nverts;++n){
      mesh->EToV[e*mesh->Nverts + n] = acceptedPartition[e].v[n];
      mesh->EX  [e*mesh->Nverts + n] = acceptedPartition[e].EX[n];
      mesh->EY  [e*mesh->Nverts + n] = acceptedPartition[e].EY[n];
      mesh->EZ  [e*mesh->Nverts + n] = acceptedPartition[e].EZ[n];
    }
  }

  // connect elements using parallel sort
  meshParallelConnect(mesh);

  // print out connectivity statistics
  meshPartitionStatistics(mesh);
  
  // connect elements to boundary faces
  meshConnectBoundary(mesh);

  // compute physical (x,y) locations of the element nodes
  meshPhysicalNodesTet3D(mesh);

  // compute geometric factors
  meshGeometricFactorsTet3D(mesh);

  // set up halo exchange info for MPI (do before connect face nodes)
  meshHaloSetup(mesh);

  // connect face nodes (find trace indices)
  meshConnectFaceNodes3D(mesh);

  // compute surface geofacs
  meshSurfaceGeometricFactorsTet3D(mesh);

  // global nodes
  meshParallelConnectNodes(mesh);

  mesh->MRABlevel = (iint *) realloc(mesh->MRABlevel,mesh->Nelements+mesh->totalHaloPairs*sizeof(iint));
  for(iint e=0;e<mesh->Nelements;++e) 
    mesh->MRABlevel[e] = acceptedPartition[e].level;
  
  free(acceptedPartition);
}

int compareCluster(const void *a, const void *b) {
  cElement_t *na = (cElement_t *) a;
  cElement_t *nb = (cElement_t *) b;

  if (na->cRank < nb->cRank) return -1;
  if (nb->cRank < na->cRank) return +1;

  if (na->cId < nb->cId) return -1;
  if (nb->cId < na->cId) return +1;  

  return 0;
}


void meshBuildMRABClusters(mesh_t *mesh, iint lev, dfloat *weights, iint *levels,
            iint *Nclusters, cluster_t **clusters, iint *Nelements, cElement_t **elements) {

  iint rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // minimum {vertex id % size}
  iint *Nsend = (iint*) calloc(size, sizeof(iint));
  iint *Nrecv = (iint*) calloc(size, sizeof(iint));
  iint *Ncount = (iint*) calloc(size, sizeof(iint));
  iint *sendOffsets = (iint*) calloc(size, sizeof(iint));
  iint *recvOffsets = (iint*) calloc(size, sizeof(iint));
  iint *sendCounts = (iint*) calloc(size, sizeof(iint));


  //build element struct
  *elements = (cElement_t *) calloc(mesh->Nelements+mesh->totalHaloPairs,sizeof(cElement_t));
  for (iint e=0;e<mesh->Nelements;e++) {
    (*elements)[e].id = e;
    (*elements)[e].level = 0.;
    if (levels) (*elements)[e].level = levels[e];
    
    (*elements)[e].weight = 1.;
    if (weights) (*elements)[e].weight = weights[e];

    for(iint n=0;n<mesh->Nverts;++n){
      (*elements)[e].v[n] = mesh->EToV[e*mesh->Nverts+n];
      (*elements)[e].EX[n] = mesh->EX[e*mesh->Nverts+n];
      (*elements)[e].EY[n] = mesh->EY[e*mesh->Nverts+n];
      (*elements)[e].EZ[n] = mesh->EZ[e*mesh->Nverts+n];
    }

    //initialize the clustering numbering
    (*elements)[e].cId = e; 
    (*elements)[e].cRank = rank;
  }

  cElement_t *sendBuffer = (cElement_t *) calloc(mesh->totalHaloPairs,sizeof(cElement_t));

  //propagate clusters 
  int allDone = 0;
  int rankDone, done;
  while(!allDone) {
    meshHaloExchange(mesh, sizeof(cElement_t), *elements, sendBuffer, *elements + mesh->Nelements);

    rankDone = 1;
    //local clustering
    done = 0;
    while(!done) {
      done = 1;
      for (iint e=0;e<mesh->Nelements;e++) {
        for (iint f=0;f<mesh->Nfaces;f++) {
          iint eP = mesh->EToE[e*mesh->Nfaces +f];
          if (eP>-1) {
            if (((*elements)[eP].level<lev+1)||((*elements)[e].level<lev+1)){
              if (compareCluster(*elements+eP,*elements+e)<0) {
                (*elements)[e].cRank = (*elements)[eP].cRank;
                (*elements)[e].cId   = (*elements)[eP].cId;
                done = 0;
                rankDone = 0;
              }
            }
          }
        }
      }    
    }

    MPI_Allreduce(&rankDone, &allDone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    allDone /= size;
  }

  //clusters have been built
  //transfer them to their owning rank

  qsort((*elements), mesh->Nelements, sizeof(cElement_t), compareCluster);

  //set up exchange along MPI interfaces
  for (iint r=0;r<size;r++)
    Nsend[r] = 0;

  for(iint e=0;e<mesh->Nelements;++e)
    ++Nsend[(*elements)[e].cRank];

  // find send offsets
  sendOffsets[0] = 0;
  for(iint r=1;r<size;++r)
    sendOffsets[r] = sendOffsets[r-1] + Nsend[r-1];
  
  // exchange byte counts 
  MPI_Alltoall(Nsend, 1, MPI_IINT,
         Nrecv, 1, MPI_IINT,
         MPI_COMM_WORLD);
  
  // count incoming faces
  iint allNrecv = 0;
  for(iint r=0;r<size;++r){
    allNrecv += Nrecv[r];
    Nrecv[r] *= sizeof(cElement_t);
    Nsend[r] *= sizeof(cElement_t);
    sendOffsets[r] *= sizeof(cElement_t);
  }
  for(iint r=1;r<size;++r)
    recvOffsets[r] = recvOffsets[r-1] + Nrecv[r-1];

  // buffer for recvied elements
  cElement_t *recvElements = (cElement_t*) calloc(allNrecv, sizeof(cElement_t));  

  // exchange parallel faces
  MPI_Alltoallv(*elements, Nsend, sendOffsets, MPI_CHAR,
      recvElements, Nrecv, recvOffsets, MPI_CHAR,
      MPI_COMM_WORLD);

  free(*elements);
  *elements = recvElements;
  *Nelements = allNrecv;

  qsort((*elements), *Nelements, sizeof(cElement_t), compareCluster);

  //build cluster lists
  // the lists are already sorted by cluster, so we just scan for different indices
  *Nclusters = 0;
  if (*Nelements) {
    (*Nclusters)++;
    for (iint e=1;e<*Nelements;e++) {
      if ((*elements)[e].cId != (*elements)[e-1].cId) (*Nclusters)++;
    }

    *clusters = (cluster_t *) calloc(*Nclusters,sizeof(cluster_t));

    iint cnt  = 0;
    iint ecnt = 1;
    (*clusters)[0].Nelements = 1;
    (*clusters)[0].offSet = 0;
    for (iint e=1;e<*Nelements;e++) {
      if ((*elements)[e].cId != (*elements)[e-1].cId) {
        cnt++;
        (*clusters)[cnt].offSet = e;
        (*clusters)[cnt].Nelements = 1;
      } else {
        (*clusters)[cnt].Nelements++;
      }
    }
  }
}

// compare the Morton indices for two clusters
int compareIndex(const void *a, const void *b){

  parallelCluster_t *ca = (parallelCluster_t*) a;
  parallelCluster_t *cb = (parallelCluster_t*) b;
  
  if(ca->index < cb->index) return -1;
  if(ca->index > cb->index) return  1;
  
  return 0;
}

// compare the Morton indices for two clusters
int compareRank(const void *a, const void *b){

  parallelCluster_t *ca = (parallelCluster_t*) a;
  parallelCluster_t *cb = (parallelCluster_t*) b;
  
  if(ca->rank < cb->rank) return -1;
  if(ca->rank > cb->rank) return  1;

  if(ca->offSet < cb->offSet) return -1;
  if(ca->offSet > cb->offSet) return  1;
  
  return 0;
}

// geometric partition of clusters of elements in 2D mesh using Morton ordering + parallelSort
dfloat meshClusteredGeometricPartition3D(mesh3D *mesh, iint Nclusters, cluster_t *clusters, 
                              iint *Nelements, cElement_t **elements){

  iint rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  iint maxNclusters;
  MPI_Allreduce(&Nclusters, &maxNclusters, 1, MPI_IINT, MPI_MAX,MPI_COMM_WORLD);
  maxNclusters = 2*((maxNclusters+1)/2);
  
  // fix maxNclusters
  parallelCluster_t *parallelClusters 
    = (parallelCluster_t*) calloc(maxNclusters, sizeof(parallelCluster_t));

  // local bounding box of element centers
  dfloat mincx = 1e9, maxcx = -1e9;
  dfloat mincy = 1e9, maxcy = -1e9;
  dfloat mincz = 1e9, maxcz = -1e9;

  // compute cluster centers on this process
  for(iint cnt=0;cnt<Nclusters;++cnt){
    iint id = clusters[cnt].offSet;
    dfloat cx = 0, cy = 0, cz = 0;
    for (iint e=0;e<clusters[cnt].Nelements;e++) {
      for(iint n=0;n<mesh->Nverts;++n){
        cx += (*elements)[id+e].EX[n];
        cy += (*elements)[id+e].EY[n];
        cz += (*elements)[id+e].EZ[n];
      }
    }
    cx /= (mesh->Nverts*clusters[cnt].Nelements);
    cy /= (mesh->Nverts*clusters[cnt].Nelements);
    cz /= (mesh->Nverts*clusters[cnt].Nelements);
    
    mincx = mymin(mincx, cx);
    maxcx = mymax(maxcx, cx);
    mincy = mymin(mincy, cy);
    maxcy = mymax(maxcy, cy);
    mincz = mymin(mincz, cz);
    maxcz = mymax(maxcz, cz);
  }

  dfloat delta = 1e-4;
  mincx -= delta;
  mincy -= delta;
  mincz -= delta;
  maxcx += delta;
  maxcy += delta;
  maxcz += delta;

  // find global bounding box of cluster centers
  dfloat gmincx, gmincy, gmincz, gmaxcx, gmaxcy, gmaxcz;
  MPI_Allreduce(&mincx, &gmincx, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&mincy, &gmincy, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&mincz, &gmincz, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&maxcx, &gmaxcx, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&maxcy, &gmaxcy, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&maxcz, &gmaxcz, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);

  // choose sub-range of Morton lattice coordinates to embed cluster centers in
  unsigned long long int Nboxes = (((unsigned long long int)1)<<(bitRange-1));
  
  // compute Morton index for each cluster
  for(iint cnt=0;cnt<Nclusters;++cnt){
    // cluster center coordinates
    dfloat cx = 0, cy = 0, cz = 0;
    parallelClusters[cnt].weight = 0.;
    iint id = clusters[cnt].offSet;
    for (iint e=0;e<clusters[cnt].Nelements;e++) {
      for(iint n=0;n<mesh->Nverts;++n){
        cx += (*elements)[id+e].EX[n];
        cy += (*elements)[id+e].EY[n];
        cz += (*elements)[id+e].EZ[n];
      }
      parallelClusters[cnt].weight += (*elements)[id+e].weight;
    }
    cx /= (mesh->Nverts*clusters[cnt].Nelements);
    cy /= (mesh->Nverts*clusters[cnt].Nelements);
    cz /= (mesh->Nverts*clusters[cnt].Nelements);

    unsigned int ix = (cx-gmincx)*Nboxes/(gmaxcx-gmincx);
    unsigned int iy = (cy-gmincy)*Nboxes/(gmaxcy-gmincy);
    unsigned int iz = (cy-gmincy)*Nboxes/(gmaxcy-gmincy);

    //fill the parallel cluster struct
    parallelClusters[cnt].index =  mortonIndex3D(ix, iy,iz);
    parallelClusters[cnt].Nelements = clusters[cnt].Nelements;
    parallelClusters[cnt].offSet = clusters[cnt].offSet;
    parallelClusters[cnt].rank = rank;
  }

  // pad cluster array with dummy clusters
  for(iint n=Nclusters;n<maxNclusters;++n){
    parallelClusters[n].Nelements = -1;
    parallelClusters[n].index = mortonIndex3D(Nboxes+1, Nboxes+1,Nboxes+1);
  }

  // odd-even parallel sort of cluster capsules based on their Morton index
  parallelSort(maxNclusters, parallelClusters, sizeof(parallelCluster_t),
         compareIndex, bogusMatch);

  iint newNclusters =0;
  for (iint n=0;n<maxNclusters;n++)
    newNclusters += (parallelClusters[n].Nelements != -1);

  //Do an initial partitioning
  dfloat localTotalWeight = 0.;
  for (iint n=0; n<newNclusters; n++) 
    localTotalWeight += parallelClusters[n].weight;

  dfloat *totalWeights = (dfloat *) calloc(size,sizeof(dfloat));
  dfloat *weightOffsets = (dfloat *) calloc(size+1,sizeof(dfloat));
  
  MPI_Allgather(&localTotalWeight, 1, MPI_DFLOAT, totalWeights, 1, MPI_DFLOAT, MPI_COMM_WORLD);

  for (iint r=0; r<size; r++)
    weightOffsets[r+1] = weightOffsets[r] + totalWeights[r];

  dfloat globalTotalWeight = weightOffsets[size];
  dfloat chunkSize = globalTotalWeight/((dfloat)size);

  iint *Nsend = (iint *) calloc(size, sizeof(iint));
  iint *Nrecv = (iint *) calloc(size, sizeof(iint));
  iint *Ncount = (iint *) calloc(size, sizeof(iint));
  iint *sendOffsets = (iint*) calloc(size, sizeof(iint));
  iint *recvOffsets = (iint*) calloc(size, sizeof(iint));

  //determine the destination rank based on which chunk the cluster is in
  localTotalWeight = weightOffsets[rank];
  for (iint n=0; n<newNclusters; n++) {
    iint destRank = (iint) (localTotalWeight/chunkSize);
    Nsend[destRank]++; 
    localTotalWeight += parallelClusters[n].weight;
  }

  // find send offsets
  for(iint r=1;r<size;++r)
    sendOffsets[r] = sendOffsets[r-1] + Nsend[r-1];
  
  // exchange byte counts 
  MPI_Alltoall(Nsend, 1, MPI_IINT, Nrecv, 1, MPI_IINT, MPI_COMM_WORLD);
  
  // count incoming clusters
  newNclusters = 0;
  for(iint r=0;r<size;++r){
    newNclusters += Nrecv[r];
    Nrecv[r] *= sizeof(parallelCluster_t);
    Nsend[r] *= sizeof(parallelCluster_t);
    sendOffsets[r] *= sizeof(parallelCluster_t);
  }
  for(iint r=1;r<size;++r)
    recvOffsets[r] = recvOffsets[r-1] + Nrecv[r-1];

  parallelCluster_t *tmpParallelClusters = (parallelCluster_t *) calloc(newNclusters, sizeof(parallelCluster_t));
  
  // exchange parallel clusters
  MPI_Alltoallv(parallelClusters, Nsend, sendOffsets, MPI_CHAR,
                tmpParallelClusters, Nrecv, recvOffsets, MPI_CHAR, MPI_COMM_WORLD);

  if (parallelClusters) free(parallelClusters);
  parallelClusters = tmpParallelClusters;

  //improve the partitioning by exchanging elements between neighboring prcesses
  dfloat partQuality = improveClusteredPartition(&newNclusters, &parallelClusters);

  //now that we're partitioned and (hopefully) balanced, send the elements

  // count number of elements that should end up on this process
  iint newNelements = 0;
  for(iint n=0;n<newNclusters;n++)
    newNelements += parallelClusters[n].Nelements;

  //record the destination info
  if (newNclusters) {
    parallelClusters[0].destId = 0;
    parallelClusters[0].destOffset = 0;
    parallelClusters[0].destRank = rank;
  }
  for (iint n=1; n<newNclusters; n++) {
    parallelClusters[n].destId = n;
    parallelClusters[n].destOffset = parallelClusters[n-1].destOffset + parallelClusters[n-1].Nelements;
    parallelClusters[n].destRank = rank;
  }

  //sort by original rank and offset
  qsort(parallelClusters, newNclusters, sizeof(parallelCluster_t), compareRank);

  //reset counters
  for(iint r=0;r<size;++r)
    Nsend[r] =0;

  for (iint n=0;n<newNclusters;n++) 
    Nsend[parallelClusters[n].rank]++;

  for(iint r=1;r<size;++r)
    sendOffsets[r] = sendOffsets[r-1] + Nsend[r-1];

  // exchange byte counts 
  MPI_Alltoall(Nsend, 1, MPI_IINT, Nrecv, 1, MPI_IINT, MPI_COMM_WORLD);

  for(iint r=0;r<size;++r){
    Nrecv[r] *= sizeof(parallelCluster_t);
    Nsend[r] *= sizeof(parallelCluster_t);
    sendOffsets[r] *= sizeof(parallelCluster_t);
  }
  for(iint r=1;r<size;++r)
    recvOffsets[r] = recvOffsets[r-1] + Nrecv[r-1];

  parallelCluster_t *recvParallelClusters;
  if (Nclusters) 
    recvParallelClusters = (parallelCluster_t*) calloc(Nclusters, sizeof(parallelCluster_t));


  MPI_Alltoallv(parallelClusters, Nsend, sendOffsets, MPI_CHAR,
              recvParallelClusters, Nrecv, recvOffsets, MPI_CHAR, MPI_COMM_WORLD);

  //build the array of elements to send
  cElement_t *sendElements = (cElement_t *) calloc(1,sizeof(cElement_t));
  cElement_t *recvElements = (cElement_t *) calloc(1,sizeof(cElement_t));

  if (*Nelements) sendElements = (cElement_t *) calloc(*Nelements,sizeof(cElement_t));
  if (newNelements) recvElements = (cElement_t *) calloc(newNelements,sizeof(cElement_t));

  //reset send counts
  for (iint r=0; r<size; r++)
    Nsend[r] = 0;

  for (iint n=0;n<Nclusters;n++) {
    Nsend[recvParallelClusters[n].destRank] += recvParallelClusters[n].Nelements;
  }

  // find send offsets
  for(iint r=1;r<size;++r)
    sendOffsets[r] = sendOffsets[r-1] + Nsend[r-1];

  //build the array of elements to send
  for (iint n=0;n<Nclusters;n++) {
    iint destRank = recvParallelClusters[n].destRank;
    iint cnt = recvParallelClusters[n].Nelements;

    iint sendId = sendOffsets[destRank] + Ncount[destRank];
    iint id = recvParallelClusters[n].offSet;
    memcpy(sendElements+sendId, *elements+id, cnt*sizeof(cElement_t)); 
    Ncount[destRank] += cnt;
  }
  free(recvParallelClusters);

  // exchange element counts 
  MPI_Alltoall(Nsend, 1, MPI_IINT, Nrecv, 1, MPI_IINT, MPI_COMM_WORLD);
  
  for(iint r=0;r<size;++r){
    Nrecv[r] *= sizeof(cElement_t);
    Nsend[r] *= sizeof(cElement_t);
    sendOffsets[r] *= sizeof(cElement_t);
  }
  for(iint r=1;r<size;++r)
    recvOffsets[r] = recvOffsets[r-1] + Nrecv[r-1];

  MPI_Alltoallv(sendElements, Nsend, sendOffsets, MPI_CHAR,
                recvElements, Nrecv, recvOffsets, MPI_CHAR, MPI_COMM_WORLD);
  free(sendElements);

  //write the clusters in the proper order
  cluster_t *newClusters = (cluster_t *) calloc(newNclusters,sizeof(cluster_t));
  cElement_t *newElements = (cElement_t *) calloc(newNelements,sizeof(cElement_t));
  iint cnt =0;
  for (iint n=0;n<newNclusters;n++) {
    iint id = parallelClusters[n].destId;
    newClusters[id].Nelements = parallelClusters[n].Nelements;
    newClusters[id].offSet = parallelClusters[n].destOffset;
    for (iint e=0;e<parallelClusters[n].Nelements;e++) {
      memcpy(newElements + newClusters[id].offSet+e, recvElements+cnt++, sizeof(cElement_t));
    }
  }
  free(recvElements);
  free(parallelClusters);

  *Nelements = newNelements;

  if (*elements) free(*elements);
  *elements = newElements;

  return partQuality;
};



//swap clusters between neighboring processes to try and improve the partitioning
void balance(iint rankL, iint rankR, dfloat *weightL, dfloat *weightR, 
              iint *Nclusters, parallelCluster_t **parallelClusters) {
  
  iint rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  int tag = 999;
  MPI_Request recv, send;
  MPI_Status status;

  if (rank==rankL) {
    if ( *weightL > *weightR) {
      //count number of clusters to send to proc
      iint Nsend = 0;
      for (iint cnt=*Nclusters-1;cnt>-1;cnt--) {
        dfloat w = (*parallelClusters)[cnt].weight;
        if ((*weightL-w)>=(*weightR+w)) {
          //sending this cluster improves the balance
          *weightL -= w;
          *weightR += w;
          Nsend++; 
        } else if((*weightL-w) > *weightR) { 
          //sending makes the neighbor have a higher weight, but it improves the balance
          *weightL -= w;
          *weightR += w;
          Nsend++; 
          break;
        } else {
          break;
        }
      }

      MPI_Isend(&Nsend, 1, MPI_IINT,  rankR, tag, MPI_COMM_WORLD, &send);
      MPI_Wait(&send, &status);

      if (Nsend) {
        *Nclusters -= Nsend;

        MPI_Isend((*parallelClusters) + *Nclusters, Nsend*sizeof(parallelCluster_t), MPI_CHAR,  rankR, tag, MPI_COMM_WORLD, &send);
        MPI_Wait(&send, &status);
      }
    } else if ( *weightL < *weightR) {
      iint Nrecv;
      MPI_Irecv(&Nrecv, 1, MPI_IINT,  rankR, tag, MPI_COMM_WORLD, &recv);
      MPI_Wait(&recv, &status);

      if (Nrecv) {
        parallelCluster_t *newParallelClusters = (parallelCluster_t *) calloc(*Nclusters+Nrecv,sizeof(parallelCluster_t));
        memcpy(newParallelClusters,*parallelClusters,*Nclusters*sizeof(parallelCluster_t));

        MPI_Irecv(newParallelClusters+*Nclusters, Nrecv*sizeof(parallelCluster_t), MPI_CHAR,  rankR, tag, MPI_COMM_WORLD, &recv);
        MPI_Wait(&recv, &status);
        
        for (iint n=*Nclusters;n<*Nclusters+Nrecv;n++) {
          dfloat w = newParallelClusters[n].weight;
          *weightL += w;
          *weightR -= w;
        }

        *Nclusters += Nrecv;
        free(*parallelClusters);
        *parallelClusters = newParallelClusters;
      }
    }
  } else if (rank==rankR) {
    if (*weightL < *weightR) {
      //count number of clusters to send to proc
      iint Nsend = 0;
      for (iint cnt=0;cnt<*Nclusters;cnt++) {
        dfloat w = (*parallelClusters)[cnt].weight;
        if ((*weightR-w)>=(*weightL+w)) {
          //sending this cluster improves the balance
          *weightR -= w;
          *weightL += w;
          Nsend++; 
        } else if((*weightR-w) > *weightL) { 
          //sending makes the neighbor have a higher weight, but it improves the balance
          *weightR -= w;
          *weightL += w;
          Nsend++; 
          break;
        } else {
          break;
        }
      }

      MPI_Isend(&Nsend, 1, MPI_IINT,  rankL, tag, MPI_COMM_WORLD, &send);
      MPI_Wait(&send, &status);

      if (Nsend) {
        *Nclusters -= Nsend;
        parallelCluster_t *newParallelClusters = (parallelCluster_t *) calloc(*Nclusters,sizeof(parallelCluster_t));
        memcpy(newParallelClusters,(*parallelClusters) + Nsend,*Nclusters*sizeof(parallelCluster_t));

        MPI_Isend(*parallelClusters, Nsend*sizeof(parallelCluster_t), MPI_CHAR,  rankL, tag, MPI_COMM_WORLD, &send);
        MPI_Wait(&send, &status);

        free(*parallelClusters);
        *parallelClusters = newParallelClusters;
      }
    } else if (*weightL > *weightR) {
      iint Nrecv;
      MPI_Irecv(&Nrecv, 1, MPI_IINT,  rankL, tag, MPI_COMM_WORLD, &recv);
      MPI_Wait(&recv, &status);

      if (Nrecv) {
        parallelCluster_t *tmpParallelClusters = (parallelCluster_t *) calloc(Nrecv,sizeof(parallelCluster_t));

        MPI_Irecv(tmpParallelClusters, Nrecv*sizeof(parallelCluster_t), MPI_CHAR, rankL, tag, MPI_COMM_WORLD, &recv);
        MPI_Wait(&recv, &status);

        for (iint n=0;n<Nrecv;n++) {
          dfloat w = tmpParallelClusters[n].weight;
          *weightR += w;
          *weightL -= w;
        }

        *Nclusters += Nrecv;
        parallelCluster_t *newParallelClusters = (parallelCluster_t *) calloc(*Nclusters,sizeof(parallelCluster_t));
        memcpy(newParallelClusters,tmpParallelClusters,Nrecv*sizeof(parallelCluster_t));
        memcpy(newParallelClusters+Nrecv,*parallelClusters,(*Nclusters-Nrecv)*sizeof(parallelCluster_t));

        free(tmpParallelClusters);
        free(*parallelClusters);
        *parallelClusters = newParallelClusters;
      }
    }
  }
}


dfloat improveClusteredPartition(iint *Nclusters, parallelCluster_t **parallelClusters){

  iint rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int tag = 999;

  MPI_Request recv, send;
  MPI_Status status;

  dfloat *totalWeights = (dfloat *) calloc(size,sizeof(dfloat));
  dfloat quality;

  while (true) {

    dfloat localTotalWeight = 0.;
    for (iint n=0; n<*Nclusters; n++) 
      localTotalWeight += (*parallelClusters)[n].weight;
    
    MPI_Allgather(&localTotalWeight, 1, MPI_DFLOAT, totalWeights, 1, MPI_DFLOAT, MPI_COMM_WORLD);

    dfloat maxTotalWeight, minTotalWeight;
    MPI_Allreduce(&localTotalWeight, &minTotalWeight, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&localTotalWeight, &maxTotalWeight, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);

    quality = minTotalWeight/maxTotalWeight;

    //ends
    if ((rank==0)||(rank==size-1)) 
      balance(size-1,0,totalWeights+size-1, totalWeights+0, Nclusters,parallelClusters);

    //resync
    localTotalWeight = totalWeights[rank];
    MPI_Allgather(&localTotalWeight, 1, MPI_DFLOAT, totalWeights, 1, MPI_DFLOAT, MPI_COMM_WORLD);

    //evens
    if (( (rank%2) == 0)&&(rank+1 < size))
      balance(rank,rank+1,totalWeights+rank, totalWeights+rank+1, Nclusters,parallelClusters);
    if (( (rank%2) == 1)&&(rank-1 > -1))
      balance(rank-1,rank,totalWeights+rank-1, totalWeights+rank, Nclusters,parallelClusters);

    //resync
    localTotalWeight = totalWeights[rank];
    MPI_Allgather(&localTotalWeight, 1, MPI_DFLOAT, totalWeights, 1, MPI_DFLOAT, MPI_COMM_WORLD);

    //odds
    if (((rank%2) == 0)&&(rank-1 > -1))
      balance(rank-1,rank,totalWeights+rank-1, totalWeights+rank, Nclusters,parallelClusters);
    if (((rank%2) == 1)&&(rank+1 < size))
      balance(rank,rank+1,totalWeights+rank, totalWeights+rank+1, Nclusters,parallelClusters);

    //resync
    localTotalWeight = totalWeights[rank];
    MPI_Allgather(&localTotalWeight, 1, MPI_DFLOAT, totalWeights, 1, MPI_DFLOAT, MPI_COMM_WORLD);
    MPI_Allreduce(&localTotalWeight, &minTotalWeight, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&localTotalWeight, &maxTotalWeight, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);

    dfloat newQuality = minTotalWeight/maxTotalWeight;

    if (newQuality == quality) break; //no change
  }

  return quality;
} 