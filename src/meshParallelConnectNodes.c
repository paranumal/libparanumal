#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh.h"

typedef struct{

  int element; // local element id
  int node;    // local node id
  int rank;    // rank of original node
  int id;      // original id
  int tag;   // original bc tag
  int haloFlag;
  
  // info on base node (lowest rank node)
  int baseElement; 
  int baseNode;
  int baseRank;
  int baseId;

  // priority tag
  int priorityTag;   // minimum (non-zero bc tag)

}parallelNode_t;

// compare on base rank for sorting suitable for destination
int parallelCompareBaseNodes(const void *a, const void *b){

  parallelNode_t *fa = (parallelNode_t*) a;
  parallelNode_t *fb = (parallelNode_t*) b;

  if(fa->baseRank < fb->baseRank) return -1;
  if(fa->baseRank > fb->baseRank) return +1;

  if(fa->baseId < fb->baseId) return -1;
  if(fa->baseId > fb->baseId) return +1;

  return 0;

}

// iteratively find a gather numbering for all local element nodes
void meshParallelConnectNodes(mesh_t *mesh){

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int localNodeCount = mesh->Np*mesh->Nelements;
  int *allLocalNodeCounts = (int*) calloc(size, sizeof(int));

  MPI_Allgather(&localNodeCount, 1, MPI_int,
		allLocalNodeCounts, 1, MPI_int,
		MPI_COMM_WORLD);
  
  int gatherNodeStart = 0;
  for(int r=0;r<rank;++r)
    gatherNodeStart += allLocalNodeCounts[r];
  
  // form continuous node numbering (local=>virtual gather)
  parallelNode_t *gatherNumbering =
    (parallelNode_t*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,
			     sizeof(parallelNode_t));

  // use local numbering
  for(int e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      int id = e*mesh->Np+n;
      gatherNumbering[id].element = e;
      gatherNumbering[id].node = n;
      gatherNumbering[id].rank = rank;
      gatherNumbering[id].id = 1 + id + mesh->Nnodes + gatherNodeStart;

      gatherNumbering[id].baseElement = e;
      gatherNumbering[id].baseNode = n;
      gatherNumbering[id].baseRank = rank;
      gatherNumbering[id].baseId = 1 + id + mesh->Nnodes + gatherNodeStart;

    }

    // use vertex ids for vertex nodes to reduce iterations
    for(int v=0;v<mesh->Nverts;++v){
      int id = e*mesh->Np + mesh->vertexNodes[v];
      int gid = mesh->EToV[e*mesh->Nverts+v] + 1;
      gatherNumbering[id].id = gid; 
      gatherNumbering[id].baseId = gid; 
    }

    // label halo flags
    for(int f=0;f<mesh->Nfaces;++f){
      if(mesh->EToP[e*mesh->Nfaces+f]!=-1){
      	for(int n=0;n<mesh->Nfp;++n){
      	  int id = e*mesh->Np+mesh->faceNodes[f*mesh->Nfp+n];
      	  gatherNumbering[id].haloFlag = 1;
      	}
      }
    }
		      
    // use element-to-boundary connectivity to create tag for local nodes
    for(int f=0;f<mesh->Nfaces;++f){
      for(int n=0;n<mesh->Nfp;++n){
      	int tag = mesh->EToB[e*mesh->Nfaces+f];
      	int id = e*mesh->Np+mesh->faceNodes[f*mesh->Nfp+n];
      	if(tag>0){
      	  gatherNumbering[id].tag = tag;
      	  gatherNumbering[id].priorityTag = tag;
      	}
      }
    }
  }

  int localChange = 0, gatherChange = 1;

  parallelNode_t *sendBuffer =
    (parallelNode_t*) calloc(mesh->totalHaloPairs*mesh->Np, sizeof(parallelNode_t));

  // keep comparing numbers on positive and negative traces until convergence
  while(gatherChange>0){

    // reset change counter
    localChange = 0;

    // send halo data and recv into extension of buffer
    meshHaloExchange(mesh, mesh->Np*sizeof(parallelNode_t),
		     gatherNumbering, sendBuffer, gatherNumbering+localNodeCount);

    // compare trace nodes
    for(int e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Nfp*mesh->Nfaces;++n){
      	int id  = e*mesh->Nfp*mesh->Nfaces + n;
      	int idM = mesh->vmapM[id];
      	int idP = mesh->vmapP[id];
      	int gidM = gatherNumbering[idM].baseId;
      	int gidP = gatherNumbering[idP].baseId;

      	int baseRankM = gatherNumbering[idM].baseRank;
      	int baseRankP = gatherNumbering[idP].baseRank;

      	// use minimum of trace variables
      	int haloM = gatherNumbering[idM].haloFlag;
      	int haloP = gatherNumbering[idP].haloFlag;
      	gatherNumbering[idM].haloFlag = mymax(haloM, haloP);
      	gatherNumbering[idP].haloFlag = mymax(haloM, haloP);
      	
      	if(gidM<gidP || (gidP==gidM && baseRankM<baseRankP)){
      	  ++localChange;
      	  gatherNumbering[idP].baseElement = gatherNumbering[idM].baseElement;
      	  gatherNumbering[idP].baseNode    = gatherNumbering[idM].baseNode;
      	  gatherNumbering[idP].baseRank    = gatherNumbering[idM].baseRank;
      	  gatherNumbering[idP].baseId      = gatherNumbering[idM].baseId;
      	}
      	
      	if(gidP<gidM || (gidP==gidM && baseRankP<baseRankM)){
      	  ++localChange;
      	  gatherNumbering[idM].baseElement = gatherNumbering[idP].baseElement;
      	  gatherNumbering[idM].baseNode    = gatherNumbering[idP].baseNode;
      	  gatherNumbering[idM].baseRank    = gatherNumbering[idP].baseRank;
      	  gatherNumbering[idM].baseId      = gatherNumbering[idP].baseId;
      	}
      	
      	int tagM = gatherNumbering[idM].priorityTag;
      	int tagP = gatherNumbering[idP].priorityTag;

      	// use maximum non-zero tag for both nodes
      	if(tagM!=tagP){
      	  if(tagP>tagM){
      	    ++localChange;
      	    gatherNumbering[idM].priorityTag = tagP;
      	  }
      	  if(tagM>tagP){
      	    ++localChange;
      	    gatherNumbering[idP].priorityTag = tagM;
      	  }
      	}
      }
    }

    // sum up changes
    MPI_Allreduce(&localChange, &gatherChange, 1, MPI_int, MPI_SUM, MPI_COMM_WORLD);

    // report
    if(rank==0)
      printf("gatherChange=%d\n", gatherChange);
  }

  // sort based on base nodes (rank,element,node at base)
  qsort(gatherNumbering, localNodeCount, sizeof(parallelNode_t), parallelCompareBaseNodes);

  // extract base index of each node (i.e. gather numbering)
  mesh->gatherLocalIds  = (int*) calloc(localNodeCount, sizeof(int));
  mesh->gatherBaseIds   = (int*) calloc(localNodeCount, sizeof(int));
  mesh->gatherBaseRanks = (int*) calloc(localNodeCount, sizeof(int));
  mesh->gatherHaloFlags = (int*) calloc(localNodeCount, sizeof(int));

  for(int id=0;id<localNodeCount;++id){
    mesh->gatherLocalIds[id]  = gatherNumbering[id].element*mesh->Np+gatherNumbering[id].node;
    mesh->gatherBaseIds[id]   = gatherNumbering[id].baseId;
    mesh->gatherBaseRanks[id] = gatherNumbering[id].baseRank;
    mesh->gatherHaloFlags[id] = gatherNumbering[id].haloFlag;
  }

  // also need to extract bc tag above !!!
  
  // should do something with tag and gather numbering arrays
  free(sendBuffer);
  free(gatherNumbering);
}
