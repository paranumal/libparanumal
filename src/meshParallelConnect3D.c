#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mesh3D.h"


typedef struct {
  int v[4]; // 4 is max nodes per face
  int element, face, rank;
  int elementN, faceN, rankN; // N for neighbor
  int NfaceVertices;
}parallelFace_t;



// comparison function that orders vertices 
// based on their combined vertex indices
int parallelCompareVertices(const void *a, 
			    const void *b){

  parallelFace_t *fa = (parallelFace_t*) a;
  parallelFace_t *fb = (parallelFace_t*) b;

  for(int n=0;n<fa->NfaceVertices;++n){
    if(fa->v[n] < fb->v[n]) return -1;
    if(fa->v[n] > fb->v[n]) return +1;
  }

  return 0;

}
/* comparison function that orders element/face
   based on their indexes */
int parallelCompareFaces(const void *a, 
			 const void *b){

  parallelFace_t *fa = (parallelFace_t*) a;
  parallelFace_t *fb = (parallelFace_t*) b;

  if(fa->rank < fb->rank) return -1;
  if(fa->rank > fb->rank) return +1;

  if(fa->element < fb->element) return -1;
  if(fa->element > fb->element) return +1;

  if(fa->face < fb->face) return -1;
  if(fa->face > fb->face) return +1;
  
  return 0;

}

void parallelFaceMatch(void *a, void *b){
  
  parallelFace_t *pfa = (parallelFace_t*) a;
  parallelFace_t *pfb = (parallelFace_t*) b;

  pfa->elementN = pfb->element;
  pfa->faceN = pfb->face;
  pfa->rankN = pfb->rank;

  pfb->elementN = pfa->element;
  pfb->faceN = pfa->face;
  pfb->rankN = pfa->rank;
}

void dummyMatch(void *, void *){ }

// mesh is the local partition
void meshParallelConnect3D(mesh3D *mesh){

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // serial connectivity on each process
  meshConnect(mesh);

  // count # of broken links
  int Nbroken =0;
  for(int e=0;e<mesh->Nelements;++e)
    for(int f=0;f<mesh->Nfaces;++f)
      if(mesh->EToE[e*mesh->Nfaces+f]==-1)
	++Nbroken;
  
  Nbroken = 2*((Nbroken+1)/2);
  
  int maxNbroken;
  MPI_Allreduce(&Nbroken,&maxNbroken, 1, MPI_int, 
		MPI_MAX, MPI_COMM_WORLD);
  
  //
  parallelFace_t *brokenFaces = 
    (parallelFace_t*) calloc(maxNbroken, sizeof(parallelFace_t));
  
  int cnt = 0;
  for(int e=0;e<mesh->Nelements;++e)
    for(int f=0;f<mesh->Nfaces;++f)
      if(mesh->EToE[e*mesh->Nfaces+f]==-1){

	// use faceVertices
	for(int n=0;n<mesh->NfaceVertices;++n)
	  brokenFaces[cnt].v[n] =
	    mesh->EToV[e*mesh->Nverts + mesh->faceVertices[f*mesh->NfaceVertices+n]];
	
	mysort(brokenFaces[cnt].v,mesh->NfaceVertices, "descending");

	brokenFaces[cnt].NfaceVertices = mesh->NfaceVertices;
	
	brokenFaces[cnt].element = e;
	brokenFaces[cnt].face = f;
	brokenFaces[cnt].rank = rank;

	brokenFaces[cnt].elementN = -1;
	brokenFaces[cnt].faceN = -1;
	brokenFaces[cnt].rankN = -1;

	++cnt;
      }

  while(cnt<maxNbroken){

    brokenFaces[cnt].NfaceVertices = mesh->NfaceVertices; // could use this

    for(int n=0;n<mesh->NfaceVertices;++n)
      brokenFaces[cnt].v[n] = -1-n;
    
    brokenFaces[cnt].element = -1;
    brokenFaces[cnt].face = -1;
    brokenFaces[cnt].rank = rank;
    
    brokenFaces[cnt].elementN = -1;
    brokenFaces[cnt].faceN = -1;
    brokenFaces[cnt].rankN = -1;
    ++cnt;
  }

  // fixed this - needed to sort maxNbroken not Nbroken
  parallelSort(maxNbroken, brokenFaces, sizeof(parallelFace_t),
	       parallelCompareVertices, parallelFaceMatch);

  parallelSort(maxNbroken, brokenFaces, sizeof(parallelFace_t),
	       parallelCompareFaces, dummyMatch);

  // extract connectivity info
  mesh->EToP = (int*) calloc(mesh->Nelements*mesh->Nfaces, sizeof(int));
  for(cnt=0;cnt<mesh->Nelements*mesh->Nfaces;++cnt)
    mesh->EToP[cnt] = -1;

  for(cnt=0;cnt<maxNbroken;++cnt){
    int e = brokenFaces[cnt].element;
    int f = brokenFaces[cnt].face;
    int eN = brokenFaces[cnt].elementN;
    int fN = brokenFaces[cnt].faceN;
    int rN = brokenFaces[cnt].rankN;

    if(e>=0 && f>=0 && eN>=0 && fN>=0){
      mesh->EToE[e*mesh->Nfaces+f] = eN;
      mesh->EToF[e*mesh->Nfaces+f] = fN;
      mesh->EToP[e*mesh->Nfaces+f] = rN;
    }
  }
}
