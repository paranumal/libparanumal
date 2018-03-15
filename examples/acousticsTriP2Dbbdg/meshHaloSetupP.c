#include <stdio.h>
#include "mpi.h"

#include "mesh.h"

typedef struct {
  
  int element, face;
  int elementN, faceN, rankN;

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
  mesh->NhaloPairs = (int*) calloc(size, sizeof(int));
  for(int e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      int r = mesh->EToP[e*mesh->Nfaces+f]; // rank of neighbor
      if(r!=-1){
        mesh->totalHaloPairs += 1;
        mesh->NhaloPairs[r] += 1;
      }
    }
  }

  // count number of MPI messages in halo exchange
  mesh->NhaloMessages = 0;
  for(int r=0;r<size;++r)
    if(mesh->NhaloPairs[r])
      ++mesh->NhaloMessages;

  // create a list of element/faces with halo neighbor
  facePair_t *haloElements = 
    (facePair_t*) calloc(mesh->totalHaloPairs, sizeof(facePair_t));

  int cnt = 0;
  for(int e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      int ef = e*mesh->Nfaces+f;
      if(mesh->EToP[ef]!=-1){
        haloElements[cnt].element  = e;
        haloElements[cnt].face     = f;
        haloElements[cnt].elementN = mesh->EToE[ef];
        haloElements[cnt].faceN    = mesh->EToF[ef];
        haloElements[cnt].rankN    = mesh->EToP[ef];
        ++cnt;
      }
    }
  }
  
  // sort the face pairs in order the destination requires
  qsort(haloElements, mesh->totalHaloPairs, sizeof(facePair_t), compareHaloFaces);

  // record the outgoing order for elements
  mesh->haloElementList = (int*) calloc(mesh->totalHaloPairs, sizeof(int));
  for(int i=0;i<mesh->totalHaloPairs;++i){
    int e = haloElements[i].element;
    mesh->haloElementList[i] = e;
  }

  // reconnect elements to ghost elements
  // (ghost elements appended to end of local element list)
  cnt = mesh->Nelements;
  for(int r=0;r<size;++r){
    for(int e=0;e<mesh->Nelements;++e){
      for(int f=0;f<mesh->Nfaces;++f){
        int ef = e*mesh->Nfaces+f;
        if(mesh->EToP[ef]==r)
          mesh->EToE[ef] = cnt++;
      }
    }
  }

  free(haloElements);
}
