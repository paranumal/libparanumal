#include <stdlib.h>
#include <stdio.h>
#include "mesh2D.h"

// structure used to encode vertices that make 
// each face, the element/face indices, and
// the neighbor element/face indices (if any)
typedef struct{

  iint element;
  iint face;

  iint v1;
  iint v2;

  iint bctype;

}boundaryFace_t;

// comparison function that orders vertices 
// based on their combined vertex indices
int compareBoundaryFaces(const void *a, 
			    const void *b){

  boundaryFace_t *fa = (boundaryFace_t*) a;
  boundaryFace_t *fb = (boundaryFace_t*) b;

  if(fa->v1 < fb->v1) return -1;
  if(fa->v1 > fb->v1) return +1;

  if(fa->v2 < fb->v2) return -1;
  if(fa->v2 > fb->v2) return +1;

  return 0;

}


/* routine to find EToB (Element To Boundary)*/
void meshConnectBoundary2D(mesh2D *mesh){

  /* count number of boundary faces (i.e. not yet connected) */
  int bcnt = 0;
  for(iint e=0;e<mesh->Nelements;++e)
    for(iint f=0;f<mesh->Nfaces;++f)
      if(mesh->EToE[e*mesh->Nfaces+f]==-1) // || mesh->EToE[e*mesh->Nfaces+f]==e)
	++bcnt;

  /* build list of boundary faces */
  boundaryFace_t *boundaryFaces = (boundaryFace_t*) calloc(bcnt+mesh->NboundaryFaces,
							   sizeof(boundaryFace_t));

  bcnt = 0; // reset counter
  for(iint e=0;e<mesh->Nelements;++e){    
    for(iint f=0;f<mesh->Nfaces;++f){
      if(mesh->EToE[e*mesh->Nfaces+f]==-1) { // || mesh->EToE[e*mesh->Nfaces+f]==e){
	iint v1 = mesh->EToV[e*mesh->Nverts + f];
	iint v2 = mesh->EToV[e*mesh->Nverts + (f+1)%mesh->Nverts];
	
	boundaryFaces[bcnt].v1 = mymax(v1,v2);
	boundaryFaces[bcnt].v2 = mymin(v1,v2);
	boundaryFaces[bcnt].element = e;
	boundaryFaces[bcnt].face = f;
	boundaryFaces[bcnt].bctype = -1;
	++bcnt;
      }
    }
  }
  
  /* add boundary info */
  for(iint b=0;b<mesh->NboundaryFaces;++b){
    iint v1 = mesh->boundaryInfo[b*3+1];
    iint v2 = mesh->boundaryInfo[b*3+2];
    boundaryFaces[bcnt].v1 = mymax(v1,v2);
    boundaryFaces[bcnt].v2 = mymin(v1,v2);
    boundaryFaces[bcnt].element = -1;
    boundaryFaces[bcnt].face = -1;
    boundaryFaces[bcnt].bctype = mesh->boundaryInfo[b*3];
    //    printf("boundaryFaces[bcnt].bctype=%d\n", boundaryFaces[bcnt].bctype);
    ++bcnt;
  }

  /* sort boundaryFaces by their vertex number pairs */
  qsort(boundaryFaces, bcnt, sizeof(boundaryFace_t), compareBoundaryFaces);

  /* scan through sorted face lists looking for element-boundary matches */
  mesh->EToB = (iint*) calloc(mesh->Nelements*mesh->Nfaces, sizeof(iint));
  for(int n=0;n<mesh->Nelements*mesh->Nfaces;++n) mesh->EToB[n] = -1;
#if 0
  for(iint cnt=0;cnt<bcnt;++cnt){
    printf("%d [%d,%d] (%d,%d) {%d}\n",
	   cnt,
	   boundaryFaces[cnt].v1,
	   boundaryFaces[cnt].v2,
	   boundaryFaces[cnt].element,
	   boundaryFaces[cnt].face,
	   boundaryFaces[cnt].bctype);
  }
#endif

  iint matches = 0;
  for(iint cnt=0;cnt<bcnt-1;++cnt){
    if(!compareBoundaryFaces(boundaryFaces+cnt, boundaryFaces+cnt+1)){
      iint e = mymax(boundaryFaces[cnt].element, boundaryFaces[cnt+1].element);
      iint f = mymax(boundaryFaces[cnt].face,    boundaryFaces[cnt+1].face);
      mesh->EToB[e*mesh->Nfaces+f] =
	mymax(boundaryFaces[cnt].bctype, boundaryFaces[cnt+1].bctype);
      //      printf("%d,%d=>%d\n", e,f,mesh->EToB[e*mesh->Nfaces+f]);
      ++matches;
    }
  }
}

