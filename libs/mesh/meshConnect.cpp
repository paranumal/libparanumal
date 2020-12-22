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

  hlong v[4];

}face_t;


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

      std::sort(faces[cnt].v, faces[cnt].v+NfaceVertices,
                std::less<hlong>());

      faces[cnt].element = e;
      faces[cnt].face = f;

      faces[cnt].elementNeighbor= -1;
      faces[cnt].faceNeighbor = -1;

      ++cnt;
    }
  }

  /* sort faces by their vertex number pairs */
  std::sort(faces, faces+Nelements*Nfaces,
            [&](const face_t& a, const face_t& b) {
              return std::lexicographical_compare(a.v, a.v+NfaceVertices,
                                                  b.v, b.v+NfaceVertices);
            });

  /* scan through sorted face lists looking for adjacent
     faces that have the same vertex ids */
  for(cnt=0;cnt<Nelements*Nfaces-1;++cnt){
    if(std::equal(faces[cnt].v, faces[cnt].v+NfaceVertices,
                  faces[cnt+1].v)){
      // match
      faces[cnt].elementNeighbor = faces[cnt+1].element;
      faces[cnt].faceNeighbor = faces[cnt+1].face;

      faces[cnt+1].elementNeighbor = faces[cnt].element;
      faces[cnt+1].faceNeighbor = faces[cnt].face;
    }
  }

  /* resort faces back to the original element/face ordering */
  std::sort(faces, faces+Nelements*Nfaces,
            [](const face_t& a, const face_t& b) {
              if(a.element < b.element) return true;
              if(a.element > b.element) return false;

              return (a.face < b.face);
            });

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
