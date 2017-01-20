#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh3D.h"
#include "mpi.h"

extern "C"
{
  void *gsParallelGatherScatterSetup(int NuniqueBases,
				     int *gatherGlobalNodes);
  
  void gsParallelGatherScatter(void *gsh, void *v, const char *type);
}

typedef struct{

  iint element; // local element id
  iint node;    // local node id
  iint rank;    // rank of original node
  iint id;      // original id
  iint tag;   // original bc tag
  
  // info on base node (lowest rank node)
  iint baseElement; 
  iint baseNode;
  iint baseRank;
  iint baseId;
  iint maxRank;
  
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


// iteratively find a global numbering for all local element nodes
void meshParallelConnectNodesHex3D(mesh3D *mesh){

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  iint localNodeCount = mesh->Np*mesh->Nelements;
  iint *allLocalNodeCounts = (iint*) calloc(size, sizeof(iint));

  MPI_Allgather(&localNodeCount, 1, MPI_IINT,
		allLocalNodeCounts, 1, MPI_IINT,
		MPI_COMM_WORLD);
  
  iint globalNodeStart = 0;
  for(iint r=0;r<rank;++r)
    globalNodeStart += allLocalNodeCounts[r];
  
  // form continuous node numbering (local=>virtual global)
  parallelNode_t *globalNumbering =
    (parallelNode_t*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np,
			     sizeof(parallelNode_t));

  // use local numbering
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      iint id = e*mesh->Np+n;
      globalNumbering[id].element = e;
      globalNumbering[id].node = n;
      globalNumbering[id].rank = rank;
      globalNumbering[id].id = 1 + id + mesh->Nnodes + globalNodeStart;

      globalNumbering[id].baseElement = e;
      globalNumbering[id].baseNode = n;
      globalNumbering[id].baseRank = rank;
      globalNumbering[id].baseId = 1 + id + mesh->Nnodes + globalNodeStart;

      globalNumbering[id].maxRank = rank;
    }

    unsigned int hash(const unsigned int value);
    
    // use vertex ids for vertex nodes to reduce iterations
    for(iint v=0;v<mesh->Nverts;++v){
      iint id = e*mesh->Np + mesh->vertexNodes[v];
      iint gid = mesh->EToV[e*mesh->Nverts+v] + 1;
      globalNumbering[id].id = gid; 
      globalNumbering[id].baseId = gid; 
    }

    // use element-to-boundary connectivity to create tag for local nodes
    for(iint f=0;f<mesh->Nfaces;++f){
      for(iint n=0;n<mesh->Nfp;++n){
	iint tag = mesh->EToB[e*mesh->Nfaces+f];
	iint id = e*mesh->Np+mesh->faceNodes[f*mesh->Nfp+n];
	if(tag>0){
	  globalNumbering[id].tag = tag;
	  globalNumbering[id].priorityTag = tag;
	}
      }
    }
  }

  iint localChange = 0, globalChange = 1;

  parallelNode_t *sendBuffer =
    (parallelNode_t*) calloc(mesh->totalHaloPairs*mesh->Np, sizeof(parallelNode_t));

  // keep comparing numbers on positive and negative traces until convergence
  while(globalChange>0){

    // reset change counter
    localChange = 0;

    // send halo data and recv into extension of buffer
    meshHaloExchange3D(mesh, mesh->Np*sizeof(parallelNode_t),
		       globalNumbering, sendBuffer, globalNumbering+localNodeCount);

    // compare trace nodes
    for(iint e=0;e<mesh->Nelements;++e){
      for(iint n=0;n<mesh->Nfp*mesh->Nfaces;++n){
	iint id  = e*mesh->Nfp*mesh->Nfaces + n;
	iint idM = mesh->vmapM[id];
	iint idP = mesh->vmapP[id];
	iint gidM = globalNumbering[idM].baseId;
	iint gidP = globalNumbering[idP].baseId;

	iint baseRankM = globalNumbering[idM].baseRank;
	iint baseRankP = globalNumbering[idP].baseRank;

	iint maxRankM = globalNumbering[idM].maxRank;
	iint maxRankP = globalNumbering[idP].maxRank;

	globalNumbering[idM].maxRank = mymax(maxRankM, maxRankP);
	globalNumbering[idP].maxRank = mymax(maxRankM, maxRankP);
	
	// use minimum of trace variables
	if(gidM<gidP || (gidP==gidM && baseRankM<baseRankP)){
	  ++localChange;
	  globalNumbering[idP].baseElement = globalNumbering[idM].baseElement;
	  globalNumbering[idP].baseNode    = globalNumbering[idM].baseNode;
	  globalNumbering[idP].baseRank    = globalNumbering[idM].baseRank;
	  globalNumbering[idP].baseId      = globalNumbering[idM].baseId;
	}

	if(gidP<gidM || (gidP==gidM && baseRankP<baseRankM)){
	  ++localChange;
	  globalNumbering[idM].baseElement = globalNumbering[idP].baseElement;
	  globalNumbering[idM].baseNode    = globalNumbering[idP].baseNode;
	  globalNumbering[idM].baseRank    = globalNumbering[idP].baseRank;
	  globalNumbering[idM].baseId      = globalNumbering[idP].baseId;
	}

	iint tagM = globalNumbering[idM].priorityTag;
	iint tagP = globalNumbering[idP].priorityTag;

	// use maximum non-zero tag for both nodes
	if(tagM!=tagP){
	  if(tagP>tagM){
	    ++localChange;
	    globalNumbering[idM].priorityTag = tagP;
	  }
	  if(tagM>tagP){
	    ++localChange;
	    globalNumbering[idP].priorityTag = tagM;
	  }
	}
      }
    }

    // sum up changes
    MPI_Allreduce(&localChange, &globalChange, 1, MPI_IINT, MPI_SUM, MPI_COMM_WORLD);

    // report
    if(rank==0)
      printf("globalChange=%d\n", globalChange);
  }  
  
  // sort based on base nodes (rank,element,node at base)
  qsort(globalNumbering, localNodeCount, sizeof(parallelNode_t), parallelCompareBaseNodes);

  // local numbers of sorted nodes
  iint *gatherLocalNodes = (iint*) calloc(localNodeCount, sizeof(iint));
  for(iint n=0;n<localNodeCount;++n){
    gatherLocalNodes[n] = globalNumbering[n].node + mesh->Np*globalNumbering[n].element; // scrape shuffled list of local indices
  }
  
  // 1. count number of unique base nodes on this process
  iint NuniqueBases = 0; // assumes at least one base node
  for(iint n=0;n<localNodeCount;++n){
    iint test = (n==0) ? 1: (globalNumbering[n].baseId != globalNumbering[n-1].baseId);
    NuniqueBases += test;
  }
  mesh->NuniqueBases = NuniqueBases;
  
  iint *gatherGlobalNodes = (iint*) calloc(NuniqueBases, sizeof(iint)); // global labels for unique nodes
  iint *gatherNodeOffsets = (iint*) calloc(NuniqueBases+1, sizeof(iint)); // offset into sorted list of nodes

  // only finds bases
  NuniqueBases = 0; // reset counter
  for(iint n=0;n<localNodeCount;++n){
    iint test = (n==0) ? 1: (globalNumbering[n].baseId != globalNumbering[n-1].baseId);
    if(test){
      gatherGlobalNodes[NuniqueBases] = globalNumbering[n].baseId; // use this for gslib gather-scatter operation
      gatherNodeOffsets[NuniqueBases++] = n;  // increment unique base counter and record index into shuffled list of ndoes

    }
  }
  gatherNodeOffsets[NuniqueBases] = localNodeCount;

  // list of nodes to extract from DEVICE gathered array
  iint NnodeHalo = 0;
  for(iint n=0;n<NuniqueBases;++n){
    int id = gatherNodeOffsets[n];
    if(globalNumbering[id].baseRank!=globalNumbering[id].maxRank){
      ++NnodeHalo;
    }
  }
  mesh->NnodeHalo = NnodeHalo;
  
  printf("NnodeHalo = %d, NuniqueBases = %d\n", NnodeHalo, NuniqueBases);

  // set up gather-scatter of halo nodes
  mesh->gsh = NULL;
  if(mesh->NnodeHalo){  
    iint *nodeHaloIds = (iint*) calloc(NnodeHalo, sizeof(iint));
    iint *nodeHaloGlobalIds = (iint*) calloc(NnodeHalo, sizeof(iint));
    NnodeHalo = 0;
    for(iint n=0;n<NuniqueBases;++n){
      int id = gatherNodeOffsets[n];
      if(globalNumbering[id].baseRank!=globalNumbering[id].maxRank){
	nodeHaloIds[NnodeHalo] = n;
	nodeHaloGlobalIds[NnodeHalo] = globalNumbering[id].baseId;
	//	printf("halo node %d base node %d global id %d\n", NnodeHalo, n, nodeHaloGlobalIds[NnodeHalo]); 
	++NnodeHalo;
      }
    }     

    // initiate gslib gather-scatter comm pattern
    mesh->gsh = gsParallelGatherScatterSetup(NnodeHalo, nodeHaloGlobalIds);
    
    // list of halo node indices
    mesh->o_nodeHaloIds = mesh->device.malloc(NnodeHalo*sizeof(iint),  nodeHaloIds);
    
    // allocate buffer for gathering on halo nodes
    mesh->subGatherTmp = (dfloat*) calloc(NnodeHalo, sizeof(dfloat));
    mesh->o_subGatherTmp  = mesh->device.malloc(NnodeHalo*sizeof(dfloat), mesh->subGatherTmp);
  }

  mesh->o_gatherTmp         = mesh->device.malloc(NuniqueBases*sizeof(dfloat));
  mesh->o_gatherNodeOffsets = mesh->device.malloc((NuniqueBases+1)*sizeof(iint), gatherNodeOffsets);
  mesh->o_gatherLocalNodes  = mesh->device.malloc(localNodeCount*sizeof(iint),   gatherLocalNodes);

  
  // should do something with tag and global numbering arrays
  free(sendBuffer);
  free(globalNumbering);
}

