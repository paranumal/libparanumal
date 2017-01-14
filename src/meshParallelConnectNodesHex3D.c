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
  dfloat x,y,z;
  
  // info on base node (lowest rank node)
  iint baseElement; 
  iint baseNode;
  iint baseRank;
  iint baseId;
  dfloat baseX, baseY, baseZ;
  
  // priority tag
  iint priorityTag;   // minimum (non-zero bc tag)

}parallelNode_t;

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
    (parallelNode_t*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np, sizeof(parallelNode_t));

  // assume ordering (brittle) (ideally search for min dist from (-1,-1)...
  // assume ordering (brittle)
  iint vnums[8];

  vnums[0] = 0;
  vnums[1] = mesh->N;
  vnums[2] = (mesh->N+1)*(mesh->N+1)-1;
  vnums[3] = mesh->N*(mesh->N+1); 

  vnums[4] = mesh->Nq*mesh->Nq*mesh->N + 0;
  vnums[5] = mesh->Nq*mesh->Nq*mesh->N + mesh->N;
  vnums[6] = mesh->Np-1;
  vnums[7] = mesh->Nq*mesh->Nq*mesh->N + mesh->N*mesh->Nq;

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

      globalNumbering[id].x = mesh->x[id];
      globalNumbering[id].y = mesh->y[id];
      globalNumbering[id].z = mesh->z[id];

      globalNumbering[id].baseX = mesh->x[id];
      globalNumbering[id].baseY = mesh->y[id];
      globalNumbering[id].baseZ = mesh->z[id];
    }

    // use vertex ids for vertex nodes to reduce iterations
    for(iint v=0;v<mesh->Nverts;++v){
      iint id = e*mesh->Np + vnums[v];
      iint gid = mesh->EToV[e*mesh->Nverts+v] + 1;
      globalNumbering[id].id = gid;
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

  parallelNode_t *sendBuffer = (parallelNode_t*) calloc(mesh->totalHaloPairs*mesh->Np, sizeof(parallelNode_t));

  printf("totalHaloPairs = %d\n", mesh->totalHaloPairs);
  
  while(globalChange>0){

    // reset change counter
    localChange = 0;

    // send halo data and recv into extension of buffer
    meshHaloExchange3D(mesh, mesh->Np*sizeof(parallelNode_t), globalNumbering, sendBuffer, globalNumbering+mesh->Np*mesh->Nelements);

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
	
	// use minimum of trace variables
	if(gidM<gidP || (gidP==gidM && baseRankM<baseRankP)){
	  ++localChange;
	  globalNumbering[idP].baseElement = globalNumbering[idM].baseElement;
	  globalNumbering[idP].baseNode    = globalNumbering[idM].baseNode;
	  globalNumbering[idP].baseRank    = globalNumbering[idM].baseRank;
	  globalNumbering[idP].baseId      = globalNumbering[idM].baseId;

	  globalNumbering[idP].baseX    = globalNumbering[idM].baseX;
	  globalNumbering[idP].baseY    = globalNumbering[idM].baseY;
	  globalNumbering[idP].baseZ    = globalNumbering[idM].baseZ;
	}

	if(gidP<gidM || (gidP==gidM && baseRankP<baseRankM)){
	  ++localChange;
	  globalNumbering[idM].baseElement = globalNumbering[idP].baseElement;
	  globalNumbering[idM].baseNode    = globalNumbering[idP].baseNode;
	  globalNumbering[idM].baseRank    = globalNumbering[idP].baseRank;
	  globalNumbering[idM].baseId      = globalNumbering[idP].baseId;

	  globalNumbering[idM].baseX    = globalNumbering[idP].baseX;
	  globalNumbering[idM].baseY    = globalNumbering[idP].baseY;
	  globalNumbering[idM].baseZ    = globalNumbering[idP].baseZ;
	}

	iint tagM = globalNumbering[idM].priorityTag;
	iint tagP = globalNumbering[idP].priorityTag;

	// use minimum non-zer tag for both nodes
	if(tagM!=tagP){
	  if(tagM>0 || tagP>0){
	    if(tagP>tagM && tagM>0){
	      ++localChange;
	      tagP = tagM;
	    }
	    if(tagM>tagP && tagP>0){
	      ++localChange;
	      tagM = tagP;
	    }
	    globalNumbering[idM].priorityTag = tagM;
	    globalNumbering[idP].priorityTag = tagP;
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

  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      iint id = e*mesh->Np+n;
      if(globalNumbering[id].baseRank>rank){
	printf("Node: (e,n,r,id,bc) = (%d,%d,%d,%d,%d) => (%d,%d,[%d],%d,(%g,%g,%g),%d) => (%d,%d,[%d],%d,(%g,%g,%g),%d)\n",
	       e,
	       n,
	       rank,
	       id,
	       globalNumbering[id].tag,
	       
	       globalNumbering[id].element,
	       globalNumbering[id].node,
	       globalNumbering[id].rank,
	       globalNumbering[id].id,
	       globalNumbering[id].x,
	       globalNumbering[id].y,
	       globalNumbering[id].z,
	       globalNumbering[id].tag,
	       globalNumbering[id].baseElement,
	       globalNumbering[id].baseNode,
	       globalNumbering[id].baseRank,
	       globalNumbering[id].baseId,
	       globalNumbering[id].baseX,
	       globalNumbering[id].baseY,
	       globalNumbering[id].baseZ,
	       globalNumbering[id].priorityTag);
      }
    }
  }

  // should do something with tag and global numbering arrays
  free(sendBuffer);
  free(globalNumbering);
}
