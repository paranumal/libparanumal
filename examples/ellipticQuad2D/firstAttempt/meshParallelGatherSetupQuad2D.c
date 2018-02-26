#include <stdlib.h>
#include "mpi.h"
#include "mesh2D.h"

typedef struct {
  int localId;
  int baseId;
  int baseRank;
  int recvIndex;
  int index;
}gatherNode_t;

// compare 
int compareGatherNodes(const void *a, const void *b){
  const gatherNode_t *nodea = (gatherNode_t*) a;
  const gatherNode_t *nodeb = (gatherNode_t*) b;

  if(nodea->baseRank < nodeb->baseRank) return -1;
  if(nodea->baseRank > nodeb->baseRank) return +1;

  if(nodea->baseId < nodeb->baseId) return -1;
  if(nodea->baseId > nodeb->baseId) return +1;
  
  return 0;
}

int compareBaseNodes(const void *a, const void *b){
  const gatherNode_t *nodea = (gatherNode_t*) a;
  const gatherNode_t *nodeb = (gatherNode_t*) b;

  if(nodea->baseId < nodeb->baseId) return -1;
  if(nodea->baseId > nodeb->baseId) return +1;
  
  return 0;
}


void meshParallelGatherSetupQuad2D(mesh2D *mesh){
  
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  // set up local-global number pairs
  gatherNode_t *gatherNodes =
    (gatherNode_t*) calloc(mesh->Np*mesh->Nelements, sizeof(gatherNode_t));

  int *gatherSendCounts = (int*) calloc(size, sizeof(int));
  
  // extract nodes that have base on other rank
  int cnt = 0;
  int gatherLocalCount = 0;
  int gatherSendCount = 0;
  
  for(int n=0;n<mesh->Np*mesh->Nelements;++n){
    int r = mesh->baseRanks[n];
    gatherNodes[n].localId  = n;
    gatherNodes[n].baseId   = mesh->baseIds[n];
    gatherNodes[n].baseRank = r;
    if(r!=rank){
      ++gatherSendCounts[r];
      ++gatherSendCount;
    }
    else
      ++gatherLocalCount;
  }
  
  // sort based on base node index [ i.e. recipient order ]
  qsort(gatherNodes, mesh->Np*mesh->Nelements,
	sizeof(gatherNode_t), compareGatherNodes);
  
  // extract: ordered list of local node indices for gather
  int *gatherBaseIds  = (int*) calloc(gatherSendCount, sizeof(int));
  int *gatherSendIds  = (int*) calloc(gatherSendCount, sizeof(int));
  
  cnt = 0;
  for(int n=0;n<mesh->Np*mesh->Nelements;++n){
    if(gatherNodes[n].baseId!=rank){
      // ids of nodes to send
      gatherSendIds[cnt] = gatherNodes[n].localId;
      
      // ids of destination nodes 
      gatherBaseIds[cnt] = gatherNodes[n].baseId;
      
      ++cnt;
    }
  }
  int gatherNsend = cnt;
  
  // all processes tell all processes how many nodes to receive
  int *gatherRecvCounts = (int*) calloc(size,sizeof(int));
  MPI_Alltoall(gatherSendCounts, 1, MPI_INT,
	       gatherRecvCounts, 1, MPI_INT,
	       MPI_COMM_WORLD);
  
  // form arrays for all to all (variable length)
  int *gatherRecvDispls = (int*) calloc(size,sizeof(int));
  int *gatherSendDispls = (int*) calloc(size,sizeof(int));
  for(int r=1;r<size;++r){
    gatherRecvDispls[r] = gatherRecvDispls[r-1]+gatherRecvCounts[r-1];
    gatherSendDispls[r] = gatherSendDispls[r-1]+gatherSendCounts[r-1];
  }

  int gatherNrecv = 0;
  for(int r=0;r<size;++r)
    gatherNrecv += gatherRecvCounts[r];
  
  int *gatherRecvIds = (int*) calloc(gatherNrecv, sizeof(int));
  MPI_Alltoallv(gatherBaseIds, gatherSendCounts, gatherSendDispls, MPI_INT,
		gatherRecvIds, gatherRecvCounts, gatherRecvDispls, MPI_INT,
		MPI_COMM_WORLD);

  // now sort incoming data
  gatherNode_t *incoming = (gatherNode_t*) calloc(gatherNrecv, sizeof(gatherNode_t));  
  for(int n=0;n<gatherNrecv;++n){
    incoming[n].baseId = gatherRecvIds[n];
    incoming[n].index = n;
  }
  qsort(incoming, gatherNrecv, sizeof(gatherNode_t), compareBaseNodes);

  // count degree of each base for incoming data
  int *gatherRecvStarts    = (int*) calloc(gatherNrecv+1, sizeof(int));
  int *gatherRecvBaseIds   = (int*) calloc(gatherNrecv,   sizeof(int));
  int *gatherSourceDegrees = (int*) calloc(gatherNrecv,   sizeof(int));
  int *gatherRecvSourceIds = (int*) calloc(gatherNrecv,   sizeof(int));

  // extract index relative to incoming
  cnt = 1; // at least one incoming
  gatherRecvSourceIds[0] = incoming[0].index;
  gatherRecvStarts[0] = 0;
  for(int n=1;n<gatherNrecv;++n){
    gatherRecvSourceIds[n] = incoming[n].index;
    if(incoming[n].baseId!=incoming[n-1].baseId)
      gatherRecvStarts[cnt++] = n;
  }
  gatherRecvStarts[cnt] = gatherNrecv;
  int gatherNbaseRecv = cnt;
  
  for(int n=0;n<gatherNbaseRecv;++n)
    gatherRecvBaseIds[n] = incoming[gatherRecvStarts[n]].baseId;

  // find local numbering stuff here
  
  mesh->gatherNsend = gatherNsend;
  mesh->gatherNrecv = gatherNrecv;
  mesh->gatherSendIds = gatherSendIds;
  mesh->gatherRecvIds = gatherRecvIds;
  mesh->gatherSendCounts = gatherSendCounts;
  mesh->gatherRecvCounts = gatherRecvCounts;
  mesh->gatherSendDispls = gatherSendDispls;
  mesh->gatherRecvDispls = gatherRecvDispls;

  mesh->gatherRecvStarts    = gatherRecvStarts;
  mesh->gatherRecvBaseIds   = gatherRecvBaseIds;
  mesh->gatherRecvSourceIds = gatherRecvSourceIds;

  mesh->gatherNbaseRecv = gatherNbaseRecv;
  
  free(gatherNodes);
  free(gatherBaseIds);
}
