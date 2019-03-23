/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "mesh.hpp"

// structure used to encode vertices that make
// each face, the element/face indices, and
// the neighbor element/face indices (if any)
typedef struct{

  dlong element;
  int face;

  dlong elementNeighbor; // neighbor element
  int faceNeighbor;    // neighbor face

  int NfaceVertices;

  hlong v[4];

}face_t;

// comparison function that orders vertices
// based on their combined vertex indices
static int compareVertices(const void *a,
                    const void *b){

  face_t *fa = (face_t*) a;
  face_t *fb = (face_t*) b;

  for(int n=0;n<fa->NfaceVertices;++n){
    if(fa->v[n] < fb->v[n]) return -1;
    if(fa->v[n] > fb->v[n]) return +1;
  }

  return 0;

}
/* comparison function that orders element/face
   based on their indexes */
static int compareFaces(const void *a,
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
void mesh_t::Connect(){

  /* build list of faces */
  face_t *faces =
    (face_t*) calloc(Nelements*Nfaces, sizeof(face_t));

  dlong cnt = 0;
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){

      for(int n=0;n<NfaceVertices;++n){
        dlong vid = e*Nverts + faceVertices[f*NfaceVertices+n];
        faces[cnt].v[n] = EToV[vid];
      }

      mysort(faces[cnt].v, NfaceVertices, "descending");

      faces[cnt].NfaceVertices = NfaceVertices;

      faces[cnt].element = e;
      faces[cnt].face = f;

      faces[cnt].elementNeighbor= -1;
      faces[cnt].faceNeighbor = -1;

      ++cnt;
    }
  }

  /* sort faces by their vertex number pairs */
  qsort(faces,
        Nelements*Nfaces,
        sizeof(face_t),
        compareVertices);

  /* scan through sorted face lists looking for adjacent
     faces that have the same vertex ids */
  for(cnt=0;cnt<Nelements*Nfaces-1;++cnt){
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
        Nelements*Nfaces,
        sizeof(face_t),
        compareFaces);

  /* extract the element to element and element to face connectivity */
  EToE = (dlong*) calloc(Nelements*Nfaces, sizeof(dlong));
  EToF = (int*)   calloc(Nelements*Nfaces, sizeof(int  ));

  cnt = 0;
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      EToE[cnt] = faces[cnt].elementNeighbor;
      EToF[cnt] = faces[cnt].faceNeighbor;
      //      printf("EToE(%d,%d) = %d \n", e,f, EToE[cnt]);
      ++cnt;
    }
  }

  // dlong Nbcs = 0;
  // for(dlong e=0;e<Nelements;++e)
  //   for(int f=0;f<Nfaces;++f)
  //     if(EToE[e*Nfaces+f]==-1)
  //       ++Nbcs;
  //printf("Nelements = %d, Nbcs = %d\n", Nelements, Nbcs);
}
