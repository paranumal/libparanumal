/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#ifdef GLIBCXX_PARALLEL
#include <parallel/algorithm>
using __gnu_parallel::sort;
#else
using std::sort;
#endif

namespace libp {

// structure used to encode vertices that make
// each face, the element/face indices, and
// the neighbor element/face indices (if any)
typedef struct{

  dlong element;
  int face;

  hlong v[4]; // max number of face vertices

  int bctype;

}boundaryFace_t;


/* routine to find EToB (Element To Boundary)*/
void mesh_t::ConnectBoundary(){

  /* count number of boundary faces (i.e. not yet connected) */
  hlong bcnt = 0;
  for(dlong e=0;e<Nelements;++e)
    for(int f=0;f<Nfaces;++f)
      if(EToE[e*Nfaces+f]==-1) // || EToE[e*Nfaces+f]==e)
        ++bcnt;

#if 0
  printf("Nbf = %d\n", NboundaryFaces);
  printf("Nfv = %d\n", NfaceVertices);
  printf("bcnt = %d\n", bcnt);
  printf("Nelements = %d\n", Nelements);
#endif

  /* build list of boundary faces */
  memory<boundaryFace_t> boundaryFaces(bcnt+NboundaryFaces);

  bcnt = 0; // reset counter
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      if(EToE[e*Nfaces+f]==-1) {

        for(int n=0;n<NfaceVertices;++n){
          dlong vid = e*Nverts + faceVertices[f*NfaceVertices+n];
          boundaryFaces[bcnt].v[n] = EToV[vid];
        }

        std::sort(boundaryFaces[bcnt].v, boundaryFaces[bcnt].v+NfaceVertices,
                  std::less<hlong>());

        boundaryFaces[bcnt].element = e;
        boundaryFaces[bcnt].face = f;
        boundaryFaces[bcnt].bctype = -1;
        ++bcnt;
      }
    }
  }

  /* add boundary info */
  for(hlong b=0;b<NboundaryFaces;++b){

    for(int n=0;n<NfaceVertices;++n)
      boundaryFaces[bcnt].v[n] = boundaryInfo[b*(NfaceVertices+1)+n+1];

    std::sort(boundaryFaces[bcnt].v, boundaryFaces[bcnt].v+NfaceVertices,
              std::less<hlong>());

    boundaryFaces[bcnt].element = -1;
    boundaryFaces[bcnt].face = -1;
    boundaryFaces[bcnt].bctype = boundaryInfo[b*(NfaceVertices+1)];

    ++bcnt;
  }

#if 0
  for(int b=0;b<bcnt;++b){
    printf("%d: e=%d, f=%d, bc=%d, v=",
           b,
           boundaryFaces[b].element,
           boundaryFaces[b].face,
           boundaryFaces[b].bctype);
    for(int n=0;n<NfaceVertices;++n)
      printf("%d ", boundaryFaces[b].v[n]);
    printf("\n");
  }
#endif

  /* sort boundaryFaces by their vertex number pairs */
  sort(boundaryFaces.ptr(), boundaryFaces.ptr()+bcnt,
      [&](const boundaryFace_t& a, const boundaryFace_t& b) {
        return std::lexicographical_compare(a.v, a.v+NfaceVertices,
                                            b.v, b.v+NfaceVertices);
      });

  /* scan through sorted face lists looking for element-boundary matches */
  EToB.malloc(Nelements*Nfaces, -1);

  #pragma omp parallel for
  for(hlong cnt=0;cnt<bcnt-1;++cnt){
    if(std::equal(boundaryFaces[cnt].v, boundaryFaces[cnt].v+NfaceVertices,
                  boundaryFaces[cnt+1].v)){
      dlong e = std::max(boundaryFaces[cnt].element, boundaryFaces[cnt+1].element);
      int f   = std::max(boundaryFaces[cnt].face,    boundaryFaces[cnt+1].face);

      EToB[e*Nfaces+f] =
        std::max(boundaryFaces[cnt].bctype, boundaryFaces[cnt+1].bctype);
    }
  }

  o_EToB = platform.malloc<int>(EToB);

#if 0
  int cnt = 0;
  for(int e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      printf("EToE(%d,%d) = %d \n", e,f, EToE[cnt]);
      ++cnt;
    }
  }
#endif
}

} //namespace libp
