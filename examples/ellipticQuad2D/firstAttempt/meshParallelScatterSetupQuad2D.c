#include <stdlib.h>
#include "mpi.h"
#include "mesh2D.h"

typedef struct {
  int localId;
  int baseId;
  int baseRank;
  int recvIndex;
  int index;
}scatterNode_t;

// compare 
int compareScatterNodes(const void *a, const void *b){
  const scatterNode_t *nodea = (scatterNode_t*) a;
  const scatterNode_t *nodeb = (scatterNode_t*) b;

  if(nodea->baseRank < nodeb->baseRank) return -1;
  if(nodea->baseRank > nodeb->baseRank) return +1;

  if(nodea->baseId < nodeb->baseId) return -1;
  if(nodea->baseId > nodeb->baseId) return +1;
  
  return 0;
}

int compareBaseNodes(const void *a, const void *b){
  const scatterNode_t *nodea = (scatterNode_t*) a;
  const scatterNode_t *nodeb = (scatterNode_t*) b;

  if(nodea->baseId < nodeb->baseId) return -1;
  if(nodea->baseId > nodeb->baseId) return +1;
  
  return 0;
}


void meshParallelScatterSetupQuad2D(mesh2D *mesh){
  
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  // set up local-global number pairs
  scatterNode_t *scatterNodes =
    (scatterNode_t*) calloc(mesh->Np*mesh->Nelements, sizeof(scatterNode_t));

  int *scatterRecvCounts = (int*) calloc(size, sizeof(int));
  
  // extract nodes that have base on other rank
  int cnt = 0;
  int scatterLocalCount = 0;
  int scatterRecvCount = 0;
  
  for(int n=0;n<mesh->Np*mesh->Nelements;++n){
    int r = mesh->baseRanks[n];
    scatterNodes[n].localId  = n;
    scatterNodes[n].baseId   = mesh->baseIds[n];
    scatterNodes[n].baseRank = r;
    if(r!=rank){
      ++scatterRecvCounts[r];
      ++scatterRecvCount;
    }
    else
      ++scatterLocalCount;
  }
  
  // sort based on scatterNode index
  qsort(scatterNodes, mesh->Np*mesh->Nelements,
	sizeof(scatterNode_t), compareScatterNodes);
  
  // extract: ordered list of local node indices for scatter
  int *scatterBaseIds  = (int*) calloc(scatterRecvCount, sizeof(int));
  int *scatterRecvIds  = (int*) calloc(scatterRecvCount, sizeof(int));
  
  cnt = 0;
  for(int n=0;n<mesh->Np*mesh->Nelements;++n){
    if(scatterNodes[n].baseId!=rank){
      // ids of nodes to recv
      scatterRecvIds[cnt] = scatterNodes[n].localId;
      
      // ids of destination nodes 
      scatterBaseIds[cnt] = scatterNodes[n].baseId;
      
      ++cnt;
    }
  }
  int scatterNrecv = cnt;
  
  // all processes tell all processes how many nodes to receive
  int *scatterSendCounts = (int*) calloc(size,sizeof(int));
  MPI_Alltoall(scatterRecvCounts, 1, MPI_int,
	       scatterSendCounts, 1, MPI_int,
	       MPI_COMM_WORLD);
  
  // form arrays for all to all (variable length)
  int *scatterSendDispls = (int*) calloc(size,sizeof(int));
  int *scatterRecvDispls = (int*) calloc(size,sizeof(int));
  for(int r=1;r<size;++r){
    scatterSendDispls[r] = scatterSendDispls[r-1]+scatterSendCounts[r-1];
    scatterRecvDispls[r] = scatterRecvDispls[r-1]+scatterRecvCounts[r-1];
  }

  int scatterNsend = 0;
  for(int r=0;r<size;++r)
    scatterNsend += scatterSendCounts[r];

  int *scatterSendIds = (int*) calloc(scatterNrecv, sizeof(int));
  MPI_Alltoallv(scatterBaseIds, scatterRecvCounts, scatterRecvDispls, MPI_int,
		scatterSendIds, scatterSendCounts, scatterSendDispls, MPI_int,
		MPI_COMM_WORLD);

  // now sort outgoing data
  scatterNode_t *outgoing = (scatterNode_t*) calloc(scatterNsend, sizeof(scatterNode_t));  
  for(int n=0;n<scatterNsend;++n){
    outgoing[n].baseId = scatterSendIds[n];
    outgoing[n].index = n;
  }
  qsort(outgoing, scatterNsend, sizeof(scatterNode_t), compareBaseNodes);

  // count degree of each base for incoming data
  int *scatterSendStarts    = (int*) calloc(scatterNsend+1, sizeof(int));
  int *scatterSendBaseIds   = (int*) calloc(scatterNsend,   sizeof(int));
  int *scatterSourceDegrees = (int*) calloc(scatterNsend,   sizeof(int));
  int *scatterSendSourceIds = (int*) calloc(scatterNsend,   sizeof(int));

  // extract index relative to incoming
  cnt = 1; // at least one incoming
  scatterSendSourceIds[0] = outgoing[0].index;
  scatterSendStarts[0] = 0;
  for(int n=1;n<scatterNsend;++n){
    scatterSendSourceIds[n] = outgoing[n].index;
    if(outgoing[n].baseId!=outgoing[n-1].baseId)
      scatterSendStarts[cnt++] = n;
  }
  scatterSendStarts[cnt] = scatterNsend;
  int scatterNbaseSend = cnt;
  
  for(int n=0;n<scatterNbaseSend;++n)
    scatterSendBaseIds[n] = outgoing[scatterSendStarts[n]].baseId;

  mesh->scatterNrecv = scatterNrecv;
  mesh->scatterNsend = scatterNsend;
  mesh->scatterRecvIds = scatterRecvIds;
  mesh->scatterSendIds = scatterSendIds;
  mesh->scatterSendCounts = scatterSendCounts;
  mesh->scatterRecvCounts = scatterRecvCounts;
  mesh->scatterSendDispls = scatterSendDispls;
  mesh->scatterRecvDispls = scatterRecvDispls;

  mesh->scatterSendStarts    = scatterSendStarts;
  mesh->scatterSendBaseIds   = scatterSendBaseIds;
  mesh->scatterSendSourceIds = scatterSendSourceIds;

  mesh->scatterNbaseSend = scatterNbaseSend;
  
  free(scatterNodes);
  free(scatterBaseIds);
}
