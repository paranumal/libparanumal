#include <stdlib.h>
#include "mpi.h"
#include "mesh2D.h"

typedef struct {
  iint localId;
  iint baseId;
  iint baseRank;
  iint recvIndex;
  iint index;
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

  iint *gatherSendCounts = (iint*) calloc(size, sizeof(iint));
  
  // extract nodes that have base on other rank
  iint cnt = 0;
  iint gatherLocalCount = 0;
  iint gatherSendCount = 0;
  
  for(iint n=0;n<mesh->Np*mesh->Nelements;++n){
    iint r = mesh->baseRanks[n];
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
  iint *gatherBaseIds  = (iint*) calloc(gatherSendCount, sizeof(iint));
  iint *gatherSendIds  = (iint*) calloc(gatherSendCount, sizeof(iint));
  
  cnt = 0;
  for(iint n=0;n<mesh->Np*mesh->Nelements;++n){
    if(gatherNodes[n].baseId!=rank){
      // ids of nodes to send
      gatherSendIds[cnt] = gatherNodes[n].localId;
      
      // ids of destination nodes 
      gatherBaseIds[cnt] = gatherNodes[n].baseId;
      
      ++cnt;
    }
  }
  iint gatherNsend = cnt;
  
  // all processes tell all processes how many nodes to receive
  iint *gatherRecvCounts = (iint*) calloc(size,sizeof(iint));
  MPI_Alltoall(gatherSendCounts, 1, MPI_IINT,
	       gatherRecvCounts, 1, MPI_IINT,
	       MPI_COMM_WORLD);
  
  // form arrays for all to all (variable length)
  iint *gatherRecvDispls = (iint*) calloc(size,sizeof(iint));
  iint *gatherSendDispls = (iint*) calloc(size,sizeof(iint));
  for(iint r=1;r<size;++r){
    gatherRecvDispls[r] = gatherRecvDispls[r-1]+gatherRecvCounts[r-1];
    gatherSendDispls[r] = gatherSendDispls[r-1]+gatherSendCounts[r-1];
  }

  iint gatherNrecv = 0;
  for(iint r=0;r<size;++r)
    gatherNrecv += gatherRecvCounts[r];
  
  iint *gatherRecvIds = (iint*) calloc(gatherNrecv, sizeof(iint));
  MPI_Alltoallv(gatherBaseIds, gatherSendCounts, gatherSendDispls, MPI_IINT,
		gatherRecvIds, gatherRecvCounts, gatherRecvDispls, MPI_IINT,
		MPI_COMM_WORLD);

  // now sort incoming data
  gatherNode_t *incoming = (gatherNode_t*) calloc(gatherNrecv, sizeof(gatherNode_t));  
  for(iint n=0;n<gatherNrecv;++n){
    incoming[n].baseId = gatherRecvIds[n];
    incoming[n].index = n;
  }
  qsort(incoming, gatherNrecv, sizeof(gatherNode_t), compareBaseNodes);

  // count degree of each base for incoming data
  iint *gatherRecvStarts    = (iint*) calloc(gatherNrecv+1, sizeof(iint));
  iint *gatherRecvBaseIds   = (iint*) calloc(gatherNrecv,   sizeof(iint));
  iint *gatherSourceDegrees = (iint*) calloc(gatherNrecv,   sizeof(iint));
  iint *gatherRecvSourceIds = (iint*) calloc(gatherNrecv,   sizeof(iint));

  // extract index relative to incoming
  cnt = 1; // at least one incoming
  gatherRecvSourceIds[0] = incoming[0].index;
  gatherRecvStarts[0] = 0;
  for(iint n=1;n<gatherNrecv;++n){
    gatherRecvSourceIds[n] = incoming[n].index;
    if(incoming[n].baseId!=incoming[n-1].baseId)
      gatherRecvStarts[cnt++] = n;
  }
  gatherRecvStarts[cnt] = gatherNrecv;
  iint gatherNbaseRecv = cnt;
  
  for(iint n=0;n<gatherNbaseRecv;++n)
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
