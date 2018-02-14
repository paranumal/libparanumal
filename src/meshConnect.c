#include <stdlib.h>
#include <stdio.h>
#include "mesh.h"

// structure used to encode vertices that make 
// each face, the element/face indices, and
// the neighbor element/face indices (if any)
typedef struct{

  iint element;
  iint face;

  iint elementNeighbor; // neighbor element
  iint faceNeighbor;    // neighbor face

  iint NfaceVertices;
  
  iint v[4];

}face_t;

// comparison function that orders vertices 
// based on their combined vertex indices
int compareVertices(const void *a, 
		    const void *b){

  face_t *fa = (face_t*) a;
  face_t *fb = (face_t*) b;
  
  for(iint n=0;n<fa->NfaceVertices;++n){
    if(fa->v[n] < fb->v[n]) return -1;
    if(fa->v[n] > fb->v[n]) return +1;
  }
  
  return 0;

}
/* comparison function that orders element/face
   based on their indexes */
int compareFaces(const void *a, 
		 const void *b){

  face_t *fa = (face_t*) a;
  face_t *fb = (face_t*) b;

  if(fa->element < fb->element) return -1;
  if(fa->element > fb->element) return +1;

  if(fa->face < fb->face) return -1;
  if(fa->face > fb->face) return +1;

  return 0;

}

/* routine to find EToE (Element To Element)
   and EToF (Element To Local Face) connectivity arrays */
void meshConnect(mesh_t *mesh){

  /* build list of faces */
  face_t *faces = 
    (face_t*) calloc(mesh->Nelements*mesh->Nfaces, sizeof(face_t));
  
  int cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint f=0;f<mesh->Nfaces;++f){
      
      for(iint n=0;n<mesh->NfaceVertices;++n){
	iint vid = e*mesh->Nverts + mesh->faceVertices[f*mesh->NfaceVertices+n];
	faces[cnt].v[n] = mesh->EToV[vid];
      }
      
      mysort(faces[cnt].v, mesh->NfaceVertices, "descending");

      faces[cnt].NfaceVertices = mesh->NfaceVertices;
      
      faces[cnt].element = e;
      faces[cnt].face = f;
      
      faces[cnt].elementNeighbor= -1;
      faces[cnt].faceNeighbor = -1;
      
      ++cnt;

    }
  }
  
  /* sort faces by their vertex number pairs */
  qsort(faces, 
	mesh->Nelements*mesh->Nfaces,
	sizeof(face_t),
	compareVertices);
  
  /* scan through sorted face lists looking for adjacent
     faces that have the same vertex ids */
  for(cnt=0;cnt<mesh->Nelements*mesh->Nfaces-1;++cnt){
    
    if(!compareVertices(faces+cnt, faces+cnt+1)){
      // match
      faces[cnt].elementNeighbor = faces[cnt+1].element;
      faces[cnt].faceNeighbor = faces[cnt+1].face;

      faces[cnt+1].elementNeighbor = faces[cnt].element;
      faces[cnt+1].faceNeighbor = faces[cnt].face;
    }
    
  }

  /* resort faces back to the original element/face ordering */
  qsort(faces, 
	mesh->Nelements*mesh->Nfaces,
	sizeof(face_t),
	compareFaces);

  /* extract the element to element and element to face connectivity */
  mesh->EToE = (iint*) calloc(mesh->Nelements*mesh->Nfaces, sizeof(iint));
  mesh->EToF = (iint*) calloc(mesh->Nelements*mesh->Nfaces, sizeof(iint));

  cnt = 0;
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint f=0;f<mesh->Nfaces;++f){
      mesh->EToE[cnt] = faces[cnt].elementNeighbor;
      mesh->EToF[cnt] = faces[cnt].faceNeighbor;
      //      printf("EToE(%d,%d) = %d \n", e,f, mesh->EToE[cnt]);
      
      ++cnt;
    }
  }

  iint Nbcs = 0;
  for(iint e=0;e<mesh->Nelements;++e)
    for(iint f=0;f<mesh->Nfaces;++f) {
      if(mesh->EToE[e*mesh->Nfaces+f]==-1)
	++Nbcs;
    }

  //printf("Nelements = %d, Nbcs = %d\n", mesh->Nelements, Nbcs);
}

