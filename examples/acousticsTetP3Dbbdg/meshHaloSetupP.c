#include <stdio.h>
#include "mpi.h"

#include "mesh.h"

typedef struct {
  
  iint element, face, N;
  iint elementN, faceN, rankN;

}facePair_t;

/* comparison function that orders halo element/face
   based on their indexes */
int compareHaloFaces(const void *a, 
         const void *b){
  
  facePair_t *fa = (facePair_t*) a;
  facePair_t *fb = (facePair_t*) b;
  
  if(fa->rankN < fb->rankN) return -1;
  if(fa->rankN > fb->rankN) return +1;

  if(fa->elementN < fb->elementN) return -1;
  if(fa->elementN > fb->elementN) return +1;

  if(fa->faceN < fb->faceN) return -1;
  if(fa->faceN > fb->faceN) return +1;
  
  return 0;
}

// set up halo infomation for inter-processor MPI 
// exchange of trace nodes
void meshHaloSetupP(mesh_t *mesh){

  // MPI info
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // non-blocking MPI isend/irecv requests (used in meshHaloExchange)
  mesh->haloSendRequests = calloc(size, sizeof(MPI_Request));
  mesh->haloRecvRequests = calloc(size, sizeof(MPI_Request));
  
  // count number of halo element nodes to swap
  mesh->totalHaloPairs = 0;
  mesh->NhaloPairs = (iint*) calloc(size, sizeof(iint));
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint f=0;f<mesh->Nfaces;++f){
      iint r = mesh->EToP[e*mesh->Nfaces+f]; // rank of neighbor
      if(r!=-1){
  mesh->totalHaloPairs += 1;
  mesh->NhaloPairs[r] += 1;
      }
    }
  }

  // count number of MPI messages in halo exchange
  mesh->NhaloMessages = 0;
  for(iint r=0;r<size;++r)
    if(mesh->NhaloPairs[r])
      ++mesh->NhaloMessages;

  // create a list of element/faces with halo neighbor
  facePair_t *haloElements = 
    (facePair_t*) calloc(mesh->totalHaloPairs, sizeof(facePair_t));

  iint cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint f=0;f<mesh->Nfaces;++f){
      iint ef = e*mesh->Nfaces+f;
      if(mesh->EToP[ef]!=-1){
  haloElements[cnt].element  = e;
  haloElements[cnt].face     = f;
  haloElements[cnt].elementN = mesh->EToE[ef];
  haloElements[cnt].faceN    = mesh->EToF[ef];
  haloElements[cnt].rankN    = mesh->EToP[ef];
  haloElements[cnt].N        = mesh->N[e];

  ++cnt;
      }
    }
  }
  
  // sort the face pairs in order the destination requires
  qsort(haloElements, mesh->totalHaloPairs, sizeof(facePair_t), compareHaloFaces);

  // record the outgoing order for elements
  mesh->haloElementList = (iint*) calloc(mesh->totalHaloPairs, sizeof(iint));
  for(iint i=0;i<mesh->totalHaloPairs;++i){
    iint e = haloElements[i].element;
    mesh->haloElementList[i] = e;
  }

  // reconnect elements to ghost elements
  // (ghost elements appended to end of local element list)
  cnt = mesh->Nelements;
  for(iint r=0;r<size;++r){
    for(iint e=0;e<mesh->Nelements;++e){
      for(iint f=0;f<mesh->Nfaces;++f){
  iint ef = e*mesh->Nfaces+f;
  if(mesh->EToP[ef]==r)
    mesh->EToE[ef] = cnt++;
      }
    }
  }

  // create halo extension for x,y arrays
  iint totalHaloNodes = mesh->totalHaloPairs*mesh->NpMax;
  iint localNodes     = mesh->Nelements*mesh->NpMax;

  // temporary send buffer
  dfloat *sendBuffer = (dfloat*) calloc(totalHaloNodes, sizeof(dfloat));

  // extend x,y arrays to hold coordinates of node coordinates of elements in halo
  mesh->x = (dfloat*) realloc(mesh->x, (localNodes+totalHaloNodes)*sizeof(dfloat));
  mesh->y = (dfloat*) realloc(mesh->y, (localNodes+totalHaloNodes)*sizeof(dfloat));
  if(mesh->dim==3)
    mesh->z = (dfloat*) realloc(mesh->z, (localNodes+totalHaloNodes)*sizeof(dfloat));
  
  // send halo data and recv into extended part of arrays
  meshHaloExchange(mesh, mesh->NpMax*sizeof(dfloat), mesh->x, sendBuffer, mesh->x + localNodes);
  meshHaloExchange(mesh, mesh->NpMax*sizeof(dfloat), mesh->y, sendBuffer, mesh->y + localNodes);
  if(mesh->dim==3)
    meshHaloExchange(mesh, mesh->NpMax*sizeof(dfloat), mesh->z, sendBuffer, mesh->z + localNodes);   

  // grab EX,EY,EZ from halo
  mesh->EX = (dfloat*) realloc(mesh->EX, (mesh->Nelements+mesh->totalHaloPairs)*mesh->Nverts*sizeof(dfloat));
  mesh->EY = (dfloat*) realloc(mesh->EY, (mesh->Nelements+mesh->totalHaloPairs)*mesh->Nverts*sizeof(dfloat));
  if(mesh->dim==3)
    mesh->EZ = (dfloat*) realloc(mesh->EZ, (mesh->Nelements+mesh->totalHaloPairs)*mesh->Nverts*sizeof(dfloat));

  // send halo data and recv into extended part of arrays
  meshHaloExchange(mesh, mesh->Nverts*sizeof(dfloat), mesh->EX, sendBuffer, mesh->EX + mesh->Nverts*mesh->Nelements);
  meshHaloExchange(mesh, mesh->Nverts*sizeof(dfloat), mesh->EY, sendBuffer, mesh->EY + mesh->Nverts*mesh->Nelements);
  if(mesh->dim==3)
    meshHaloExchange(mesh, mesh->Nverts*sizeof(dfloat), mesh->EZ, sendBuffer, mesh->EZ + mesh->Nverts*mesh->Nelements);
  
  iint *intSendBuffer = (iint *) calloc(totalHaloNodes,sizeof(iint));

  mesh->N = (iint*) realloc(mesh->N, (mesh->Nelements+mesh->totalHaloPairs)*sizeof(iint));
  
  // send halo data and recv into extended part of arrays
  meshHaloExchange(mesh, sizeof(iint), mesh->N, sendBuffer, mesh->N + mesh->Nelements);  

  free(haloElements);
  free(sendBuffer);
  free(intSendBuffer);
}
