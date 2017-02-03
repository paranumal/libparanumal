#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh2D.h"
#include "mpi.h"

typedef struct{

  iint element; // local element id
  iint node;    // local node id
  iint rank;    // rank of original node
  iint id;      // original id
  iint tag;   // original bc tag
  iint haloFlag;
  
  // info on base node (lowest rank node)
  iint baseElement; 
  iint baseNode;
  iint baseRank;
  iint baseId;

  // priority tag
  iint priorityTag;   // minimum (non-zero bc tag)

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
void meshParallelConnectNodesQuad2D(mesh2D *mesh){

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  iint localNodeCount = mesh->Np*mesh->Nelements;
  iint *allLocalNodeCounts = (iint*) calloc(size, sizeof(iint));

  MPI_Allgather(&localNodeCount, 1, MPI_IINT,
		allLocalNodeCounts, 1, MPI_IINT,
		MPI_COMM_WORLD);
  
  iint gatherNodeStart = 0;
  for(iint r=0;r<rank;++r)
    gatherNodeStart += allLocalNodeCounts[r];
  
  // form continuous node numbering (local=>virtual gather)
  parallelNode_t *gatherNumbering =
    (parallelNode_t*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,
			     sizeof(parallelNode_t));

  // use local numbering
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      iint id = e*mesh->Np+n;
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
    for(iint v=0;v<mesh->Nverts;++v){
      iint id = e*mesh->Np + mesh->vertexNodes[v];
      iint gid = mesh->EToV[e*mesh->Nverts+v] + 1;
      gatherNumbering[id].id = gid; 
      gatherNumbering[id].baseId = gid; 
    }

    // label halo flags
    for(iint f=0;f<mesh->Nfaces;++f){
      if(mesh->EToP[e*mesh->Nfaces+f]!=-1){
	for(iint n=0;n<mesh->Nfp;++n){
	  iint id = e*mesh->Np+mesh->faceNodes[f*mesh->Nfp+n];
	  gatherNumbering[id].haloFlag = 1;
	}
      }
    }
		      
    // use element-to-boundary connectivity to create tag for local nodes
    for(iint f=0;f<mesh->Nfaces;++f){
      for(iint n=0;n<mesh->Nfp;++n){
	iint tag = mesh->EToB[e*mesh->Nfaces+f];
	iint id = e*mesh->Np+mesh->faceNodes[f*mesh->Nfp+n];
	if(tag>0){
	  gatherNumbering[id].tag = tag;
	  gatherNumbering[id].priorityTag = tag;
	}
      }
    }
  }

  iint localChange = 0, gatherChange = 1;

  parallelNode_t *sendBuffer =
    (parallelNode_t*) calloc(mesh->totalHaloPairs*mesh->Np, sizeof(parallelNode_t));

  // keep comparing numbers on positive and negative traces until convergence
  while(gatherChange>0){

    // reset change counter
    localChange = 0;

    // send halo data and recv into extension of buffer
    meshHaloExchange2D(mesh, mesh->Np*sizeof(parallelNode_t),
		       gatherNumbering, sendBuffer, gatherNumbering+localNodeCount);

    // compare trace nodes
    for(iint e=0;e<mesh->Nelements;++e){
      for(iint n=0;n<mesh->Nfp*mesh->Nfaces;++n){
	iint id  = e*mesh->Nfp*mesh->Nfaces + n;
	iint idM = mesh->vmapM[id];
	iint idP = mesh->vmapP[id];
	iint gidM = gatherNumbering[idM].baseId;
	iint gidP = gatherNumbering[idP].baseId;

	iint baseRankM = gatherNumbering[idM].baseRank;
	iint baseRankP = gatherNumbering[idP].baseRank;

	// use minimum of trace variables
	iint haloM = gatherNumbering[idM].haloFlag;
	iint haloP = gatherNumbering[idP].haloFlag;
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
	
	iint tagM = gatherNumbering[idM].priorityTag;
	iint tagP = gatherNumbering[idP].priorityTag;

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
    MPI_Allreduce(&localChange, &gatherChange, 1, MPI_IINT, MPI_SUM, MPI_COMM_WORLD);
  }

  // sort based on base nodes (rank,element,node at base)
  qsort(gatherNumbering, localNodeCount, sizeof(parallelNode_t), parallelCompareBaseNodes);

  // extract base index of each node (i.e. gather numbering)
  mesh->gatherLocalIds  = (iint*) calloc(localNodeCount, sizeof(iint));
  mesh->gatherBaseIds   = (iint*) calloc(localNodeCount, sizeof(iint));
  mesh->gatherBaseRanks = (iint*) calloc(localNodeCount, sizeof(iint));
  mesh->gatherHaloFlags = (iint*) calloc(localNodeCount, sizeof(iint));

  for(iint id=0;id<localNodeCount;++id){
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
