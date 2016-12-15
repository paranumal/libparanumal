#include <stdlib.h>
#include "mpi.h"
#include "mesh2D.h"

typedef struct{
  iint localId;
  iint localRank;
  iint baseId;
  iint baseRank;
  iint globalId;
}nodeInfo_t;

// compare 
int compareNodes(const void *a, const void *b){
  const nodeInfo_t *nodea = (nodeInfo_t*) a;
  const nodeInfo_t *nodeb = (nodeInfo_t*) b;

  if(nodea->baseRank < nodeb->baseRank) return -1;
  if(nodea->baseRank > nodeb->baseRank) return +1;

  if(nodea->baseId < nodeb->baseId) return -1;
  if(nodea->baseId > nodeb->baseId) return +1;
  
  return 0;
}



// create a global numbering of nodes
// [  reqires (EToE,EToF) and halo information ]
void meshNumberNodes2D(mesh2D *mesh){

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // number of local and halo nodes
  iint NlocalNodes = mesh->Nelements*mesh->Np;
  iint NhaloNodes  = mesh->totalHaloPairs*mesh->Np;

  // number of local nodes on each process
  iint *allNlocalNodes = (iint*) calloc(size, sizeof(iint));
  MPI_Allgather(&NlocalNodes, 1, MPI_IINT,
		allNlocalNodes, 1, MPI_IINT,
		MPI_COMM_WORLD);

  iint *nodeOffsets = (iint*) calloc(size, sizeof(iint));
  for(iint r=1;r<size;++r){
    nodeOffsets[r] = nodeOffsets[r-1] + allNlocalNodes[r-1];
  }
  
  // build initial index array
  nodeInfo_t *nodeInfo = (nodeInfo_t*)
    calloc( (NlocalNodes+NhaloNodes), sizeof(nodeInfo_t));

  // start using accumulated count of nodes on lower rank processes
  for(iint n=0;n<mesh->Np*mesh->Nelements;++n){
    nodeInfo[n].baseId  = n;
    nodeInfo[n].localId   = n;
    nodeInfo[n].baseRank = rank;
    nodeInfo[n].localRank  = rank;
    nodeInfo[n].globalId = n + nodeOffsets[rank];
  }

  // halo send buffer for globalNumbering
  nodeInfo_t *sendBuffer = (nodeInfo_t*) calloc(NhaloNodes, sizeof(nodeInfo_t));
  iint maxChanges;
  
  do{  // keep comparing face node numbers

    // swap halos between partitions
    meshHaloExchange2D(mesh,mesh->Np*sizeof(nodeInfo_t),
		       nodeInfo,sendBuffer,nodeInfo+NlocalNodes);

    // loop over face nodes
    iint changes = 0;
    for(int n=0;n<mesh->Nelements*mesh->Nfaces*mesh->Nfp;++n){
      iint idM = mesh->vmapM[n];
      iint idP = mesh->vmapP[n];
      iint rankM = nodeInfo[idM].baseRank;
      iint rankP = nodeInfo[idP].baseRank;
      iint gnM = nodeInfo[idM].globalId;
      iint gnP = nodeInfo[idP].globalId;
      // use maximum of current globalNumbering for either trace
      if(gnM>gnP){
	nodeInfo[idP].baseId   = nodeInfo[idM].baseId;
	nodeInfo[idP].baseRank = nodeInfo[idM].baseRank;
	nodeInfo[idP].globalId = gnM;
	++changes;
      }
      else if(gnP>gnM){
	nodeInfo[idM].baseId   = nodeInfo[idP].baseId;
	nodeInfo[idM].baseRank = nodeInfo[idP].baseRank;
	nodeInfo[idM].globalId = gnP;
	++changes;
      }
    }
    
    // find maximum number of changes made by any process
    MPI_Allreduce(&changes, &maxChanges, 1, MPI_IINT, MPI_MAX, MPI_COMM_WORLD);
    if(rank==0) printf("number of changes = %d\n", maxChanges);

  }while(maxChanges>0); // until no face node numbers change

  qsort(nodeInfo, mesh->Np*mesh->Nelements, sizeof(nodeInfo_t), compareNodes);

  // ok to here ----->
  mesh->baseIds   = (iint*) calloc(mesh->Nelements*mesh->Np, sizeof(iint));
  mesh->baseRanks = (iint*) calloc(mesh->Nelements*mesh->Np, sizeof(iint));
  for(iint n=0;n<mesh->Np*mesh->Nelements;++n){
    mesh->baseIds[n] = nodeInfo[n].baseId;
    mesh->baseRanks[n] = nodeInfo[n].baseRank;
  }
}


