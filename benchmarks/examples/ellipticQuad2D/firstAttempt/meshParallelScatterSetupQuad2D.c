#include <stdlib.h>
#include "mpi.h"
#include "mesh2D.h"

typedef struct {
  iint localId;
  iint baseId;
  iint baseRank;
  iint recvIndex;
  iint index;
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

  iint *scatterRecvCounts = (iint*) calloc(size, sizeof(iint));
  
  // extract nodes that have base on other rank
  iint cnt = 0;
  iint scatterLocalCount = 0;
  iint scatterRecvCount = 0;
  
  for(iint n=0;n<mesh->Np*mesh->Nelements;++n){
    iint r = mesh->baseRanks[n];
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
  iint *scatterBaseIds  = (iint*) calloc(scatterRecvCount, sizeof(iint));
  iint *scatterRecvIds  = (iint*) calloc(scatterRecvCount, sizeof(iint));
  
  cnt = 0;
  for(iint n=0;n<mesh->Np*mesh->Nelements;++n){
    if(scatterNodes[n].baseId!=rank){
      // ids of nodes to recv
      scatterRecvIds[cnt] = scatterNodes[n].localId;
      
      // ids of destination nodes 
      scatterBaseIds[cnt] = scatterNodes[n].baseId;
      
      ++cnt;
    }
  }
  iint scatterNrecv = cnt;
  
  // all processes tell all processes how many nodes to receive
  iint *scatterSendCounts = (iint*) calloc(size,sizeof(iint));
  MPI_Alltoall(scatterRecvCounts, 1, MPI_IINT,
	       scatterSendCounts, 1, MPI_IINT,
	       MPI_COMM_WORLD);
  
  // form arrays for all to all (variable length)
  iint *scatterSendDispls = (iint*) calloc(size,sizeof(iint));
  iint *scatterRecvDispls = (iint*) calloc(size,sizeof(iint));
  for(iint r=1;r<size;++r){
    scatterSendDispls[r] = scatterSendDispls[r-1]+scatterSendCounts[r-1];
    scatterRecvDispls[r] = scatterRecvDispls[r-1]+scatterRecvCounts[r-1];
  }

  iint scatterNsend = 0;
  for(iint r=0;r<size;++r)
    scatterNsend += scatterSendCounts[r];

  iint *scatterSendIds = (iint*) calloc(scatterNrecv, sizeof(iint));
  MPI_Alltoallv(scatterBaseIds, scatterRecvCounts, scatterRecvDispls, MPI_IINT,
		scatterSendIds, scatterSendCounts, scatterSendDispls, MPI_IINT,
		MPI_COMM_WORLD);

  // now sort outgoing data
  scatterNode_t *outgoing = (scatterNode_t*) calloc(scatterNsend, sizeof(scatterNode_t));  
  for(iint n=0;n<scatterNsend;++n){
    outgoing[n].baseId = scatterSendIds[n];
    outgoing[n].index = n;
  }
  qsort(outgoing, scatterNsend, sizeof(scatterNode_t), compareBaseNodes);

  // count degree of each base for incoming data
  iint *scatterSendStarts    = (iint*) calloc(scatterNsend+1, sizeof(iint));
  iint *scatterSendBaseIds   = (iint*) calloc(scatterNsend,   sizeof(iint));
  iint *scatterSourceDegrees = (iint*) calloc(scatterNsend,   sizeof(iint));
  iint *scatterSendSourceIds = (iint*) calloc(scatterNsend,   sizeof(iint));

  // extract index relative to incoming
  cnt = 1; // at least one incoming
  scatterSendSourceIds[0] = outgoing[0].index;
  scatterSendStarts[0] = 0;
  for(iint n=1;n<scatterNsend;++n){
    scatterSendSourceIds[n] = outgoing[n].index;
    if(outgoing[n].baseId!=outgoing[n-1].baseId)
      scatterSendStarts[cnt++] = n;
  }
  scatterSendStarts[cnt] = scatterNsend;
  iint scatterNbaseSend = cnt;
  
  for(iint n=0;n<scatterNbaseSend;++n)
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
