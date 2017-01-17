#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mesh3D.h"

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

#if 0
  if(fa->baseElement < fb->baseElement) return -1;
  if(fa->baseElement > fb->baseElement) return +1;

  if(fa->baseNode < fb->baseNode) return -1;
  if(fa->baseNode > fb->baseNode) return +1;
#else
  if(fa->baseId < fb->baseId) return -1;
  if(fa->baseId > fb->baseId) return +1;
#endif
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

  // assume ordering (brittle) (ideally search for min dist from (-1,-1)...
  // assume ordering (brittle)
  // move this to meshLoadNodes [ then use generic meshParallelConnect ]
  iint vertexNodes[8];

  vertexNodes[0] = 0;
  vertexNodes[1] = mesh->N;
  vertexNodes[2] = (mesh->N+1)*(mesh->N+1)-1;
  vertexNodes[3] = mesh->N*(mesh->N+1); 

  vertexNodes[4] = mesh->Nq*mesh->Nq*mesh->N + 0;
  vertexNodes[5] = mesh->Nq*mesh->Nq*mesh->N + mesh->N;
  vertexNodes[6] = mesh->Np-1;
  vertexNodes[7] = mesh->Nq*mesh->Nq*mesh->N + mesh->N*mesh->Nq;

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

    // use vertex ids for vertex nodes to reduce iterations
    for(iint v=0;v<mesh->Nverts;++v){
      iint id = e*mesh->Np + vertexNodes[v];
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

	// use minimum non-zer tag for both nodes
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
  // at this point all nodes on this rank have bases on this rank or below

  
  // sort based on base nodes (rank,element,node at base)
  qsort(globalNumbering, localNodeCount,
	sizeof(parallelNode_t), parallelCompareBaseNodes);

  // A. gather on DEVICE [ local and non-local base nodes ]
  // B. copy non-local DEVICE>HOST [ contiguous block ]
  // C. gslib: gather-scatter between HOSTS<>HOSTS
  // D. copy back to buffer HOST>DEVICE
  // E. scatter
  
  // DANG - was sorting on rank/element/node of base nodes

  // 1. count number of unique base nodes on this process
  iint NuniqueBases = 0; // assumes at least one base node
  for(iint n=0;n<localNodeCount;++n){
    iint test = (n==0) ? 1: (globalNumbering[n].baseId != globalNumbering[n-1].baseId);
    NuniqueBases += test;
  }
  
  iint *gatherNsends = (iint*) calloc(size, sizeof(iint));

  iint *gatherGlobalNodes = (iint*) calloc(NuniqueBases, sizeof(iint)); // global labels for unique nodes
  iint *gatherNodeOffsets = (iint*) calloc(NuniqueBases+1, sizeof(iint)); // offset into sorted list of nodes
  
  // [strip out singletons in gather and stages ?]
  NuniqueBases = 0; // reset counter
  for(iint n=0;n<localNodeCount;++n){
    iint test = (n==0) ? 1: (globalNumbering[n].baseId != globalNumbering[n-1].baseId);
    if(test){
      gatherGlobalNodes[NuniqueBases] = globalNumbering[n].baseId; // use this for gslib gather-scatter operation
      gatherNodeOffsets[NuniqueBases++] = n;  // increment unique base counter and record index into shuffled list of ndoes

      ++gatherNsends[globalNumbering[n].baseRank];
    }
  }
  gatherNodeOffsets[NuniqueBases] = localNodeCount;

  iint *gatherLocalNodes = (iint*) calloc(localNodeCount, sizeof(iint));
  for(iint n=0;n<localNodeCount;++n){
    gatherLocalNodes[n] = globalNumbering[n].node + mesh->Np*globalNumbering[n].element; // scrape shuffled list of local indices
  }

  iint *gatherSendOffsets = (iint*) calloc(size+1, sizeof(iint));
  for(iint r=1;r<size+1;++r)
    gatherSendOffsets[r] = gatherSendOffsets[r-1]+gatherNsends[r-1];

  char allname[BUFSIZ], redname[BUFSIZ];

  sprintf(allname, "all_rank%04d.dat", rank);
  sprintf(redname, "red_rank%04d.dat", rank);
  
  FILE *allfp = (FILE*) fopen(allname, "w");
  FILE *redfp = (FILE*) fopen(redname, "w");

  fprintf(allfp, "all: \n (element, node, rank, id, tag) => base (element, node, rank, id), max rank, priority tag\n");
  for(iint n=0;n<localNodeCount;++n){
    parallelNode_t gn = globalNumbering[n];
    iint id = gn.element*mesh->Np+gn.node;
    fprintf(allfp,"(%f,%f,%f) (%05d %05d %02d %05d %01d) => (%05d %05d %02d id%05d), %02d, t%01d\n",
	    mesh->x[id],mesh->y[id],mesh->z[id],
	    gn.element, gn.node, gn.rank, gn.id, gn.tag,
	    gn.baseElement, gn.baseNode, gn.baseRank, gn.baseId,
	    gn.maxRank, gn.priorityTag);
  }
  fclose(allfp);
  fclose(redfp);

  // should do something with tag and global numbering arrays
  free(sendBuffer);
  free(globalNumbering);
}
