#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "mesh3D.h"
#include "mpi.h"

#if 0
extern "C"
{
  void *gsParallelGatherScatterSetup(int NuniqueBases, int *gatherGatherNodes);
  void  gsParallelGatherScatter(void *gsh, void *v, const char *type);
}
#endif

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

// iteratively find a gather numbering for all local element nodes
void meshParallelConnectNodesHex3D(mesh3D *mesh){

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

      gatherNumbering[id].maxRank = rank;
    }

    // use vertex ids for vertex nodes to reduce iterations
    for(iint v=0;v<mesh->Nverts;++v){
      iint id = e*mesh->Np + mesh->vertexNodes[v];
      iint gid = mesh->EToV[e*mesh->Nverts+v] + 1;
      gatherNumbering[id].id = gid; 
      gatherNumbering[id].baseId = gid; 
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
    meshHaloExchange3D(mesh, mesh->Np*sizeof(parallelNode_t),
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

	iint maxRankM = gatherNumbering[idM].maxRank;
	iint maxRankP = gatherNumbering[idP].maxRank;

	gatherNumbering[idM].maxRank = mymax(maxRankM, maxRankP);
	gatherNumbering[idP].maxRank = mymax(maxRankM, maxRankP);
	
	// use minimum of trace variables
	
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

    // report
    if(rank==0)
      printf("gatherChange=%d\n", gatherChange);
  }

  // sort based on base nodes (rank,element,node at base)
  qsort(gatherNumbering, localNodeCount, sizeof(parallelNode_t), parallelCompareBaseNodes);

  // extract base index of each node (i.e. gather numbering)
  mesh->gatherLocalIds  = (iint*) calloc(localNodeCount, sizeof(iint));
  mesh->gatherBaseIds   = (iint*) calloc(localNodeCount, sizeof(iint));
  mesh->gatherBaseRanks = (iint*) calloc(localNodeCount, sizeof(iint));
  mesh->gatherMaxRanks  = (iint*) calloc(localNodeCount, sizeof(iint));
  for(iint id=0;id<localNodeCount;++id){
    mesh->gatherLocalIds[id]  = gatherNumbering[id].element*mesh->Np+gatherNumbering[id].node;
    mesh->gatherBaseIds[id]   = gatherNumbering[id].baseId;
    mesh->gatherBaseRanks[id] = gatherNumbering[id].baseRank;
    mesh->gatherMaxRanks[id]  = gatherNumbering[id].maxRank;
  }

  // set up gslib MPI gather-scatter and OCCA gather/scatter arrays
  meshParallelGatherScatterSetup3D(mesh,
				   localNodeCount,
				   sizeof(dfloat),
				   mesh->gatherLocalIds,
				   mesh->gatherBaseIds, 
				   mesh->gatherBaseRanks,
				   mesh->gatherMaxRanks,
				   mesh->NuniqueBases,
				   mesh->o_gatherNodeOffsets,
				   mesh->o_gatherLocalNodes,
				   mesh->o_gatherTmp,
				   mesh->NnodeHalo,
				   mesh->o_nodeHaloIds,
				   mesh->o_subGatherTmp,
				   (void**)&(mesh->subGatherTmp),
				   (void**)&(mesh->gsh));

  
  // find maximum degree
  {
    for(iint n=0;n<mesh->Np*mesh->Nelements;++n){
      mesh->rhsq[n] = 1;
    }
    mesh->o_rhsq.copyFrom(mesh->rhsq);
    
    void meshParallelGatherScatter3D(mesh3D *mesh, occa::memory &o_v, occa::memory &o_gsv, const char *type);
    meshParallelGatherScatter3D(mesh, mesh->o_rhsq, mesh->o_rhsq, dfloatString);

    mesh->o_rhsq.copyTo(mesh->rhsq);
    
    dfloat maxDegree = 0, minDegree = 1e9;
    dfloat gatherMaxDegree = 0, gatherMinDegree = 1e9;
    for(iint n=0;n<mesh->Np*mesh->Nelements;++n){
      maxDegree = mymax(maxDegree, mesh->rhsq[n]);
      minDegree = mymin(minDegree, mesh->rhsq[n]);
    }

    MPI_Allreduce(&maxDegree, &gatherMaxDegree, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&minDegree, &gatherMinDegree, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

    if(rank==0){
      printf("max degree = " dfloatFormat "\n", gatherMaxDegree);
      printf("min degree = " dfloatFormat "\n", gatherMinDegree);
    }
  }
  
  // should do something with tag and gather numbering arrays
  free(sendBuffer);
  free(gatherNumbering);
}

