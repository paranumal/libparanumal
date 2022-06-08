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
typedef struct {
  hlong v[4]; // vertices on face
  dlong element, elementN;
  int face, faceN;    // face info
  int rank, rankN; // N for neighbor face info

}face_t;


// mesh is the local partition
void mesh_t::Connect(){

  EToE.malloc(Nelements*Nfaces);
  EToF.malloc(Nelements*Nfaces);
  EToP.malloc(Nelements*Nfaces);

  /**********************
   * Local Connectivity
   **********************/

  /* build list of faces */
  memory<face_t> faces(Nelements*Nfaces);

  #pragma omp parallel for collapse(2)
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){

      const dlong id = f + e*Nfaces;

      for(int n=0;n<NfaceVertices;++n){
        dlong vid = e*Nverts + faceVertices[f*NfaceVertices+n];
        faces[id].v[n] = EToV[vid];
      }

      std::sort(faces[id].v, faces[id].v+NfaceVertices,
                std::less<hlong>());

      faces[id].element = e;
      faces[id].face = f;

      faces[id].elementN= -1;
      faces[id].faceN = -1;
    }
  }

  /* sort faces by their vertex number pairs */
  sort(faces.ptr(), faces.ptr()+Nelements*Nfaces,
       [&](const face_t& a, const face_t& b) {
         return std::lexicographical_compare(a.v, a.v+NfaceVertices,
                                             b.v, b.v+NfaceVertices);
       });

  /* scan through sorted face lists looking for adjacent
     faces that have the same vertex ids */
  #pragma omp parallel for
  for(dlong cnt=0;cnt<Nelements*Nfaces-1;++cnt){
    if(std::equal(faces[cnt].v, faces[cnt].v+NfaceVertices,
                  faces[cnt+1].v)){
      // match
      faces[cnt].elementN = faces[cnt+1].element;
      faces[cnt].faceN = faces[cnt+1].face;

      faces[cnt+1].elementN = faces[cnt].element;
      faces[cnt+1].faceN = faces[cnt].face;
    }
  }

  /* resort faces back to the original element/face ordering */
  sort(faces.ptr(), faces.ptr()+Nelements*Nfaces,
       [](const face_t& a, const face_t& b) {
         if(a.element < b.element) return true;
         if(a.element > b.element) return false;

         return (a.face < b.face);
       });

  /* extract the element to element and element to face connectivity */
  #pragma omp parallel for collapse(2)
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      const dlong id = f + e*Nfaces;

      EToE[id] = faces[id].elementN;
      EToF[id] = faces[id].faceN;
    }
  }
  faces.free();


  /*****************************
   * Interprocess Connectivity
   *****************************/

  // count # of elements to send to each rank based on
  // minimum {vertex id % size}
  memory<int> Nsend(size, 0);
  memory<int> Nrecv(size, 0);
  memory<int> sendOffsets(size, 0);
  memory<int> recvOffsets(size, 0);

  // WARNING: In some corner cases, the number of faces to send may overrun int storage
  int allNsend = 0;
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      if(EToE[e*Nfaces+f]==-1){
        // find rank of destination for sorting based on max(face vertices)%size
        hlong maxv = 0;
        for(int n=0;n<NfaceVertices;++n){
          int nid = faceVertices[f*NfaceVertices+n];
          hlong id = EToV[e*Nverts + nid];
          maxv = std::max(maxv, id);
        }
        int destRank = (int) (maxv%size);

        // increment send size for
        ++Nsend[destRank];
        ++allNsend;
      }
    }
  }

  // find send offsets
  for(int rr=1;rr<size;++rr)
    sendOffsets[rr] = sendOffsets[rr-1] + Nsend[rr-1];

  // reset counters
  for(int rr=0;rr<size;++rr)
    Nsend[rr] = 0;

  // buffer for outgoing data
  memory<face_t> sendFaces(allNsend);

  // pack face data
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      if(EToE[e*Nfaces+f]==-1){

        // find rank of destination for sorting based on max(face vertices)%size
        hlong maxv = 0;
        for(int n=0;n<NfaceVertices;++n){
          int nid = faceVertices[f*NfaceVertices+n];
          hlong id = EToV[e*Nverts + nid];
          maxv = std::max(maxv, id);
        }
        int destRank = (int) (maxv%size);

        // populate face to send out staged in segment of sendFaces array
        int id = sendOffsets[destRank]+Nsend[destRank];

        sendFaces[id].element = e;
        sendFaces[id].face = f;
        for(int n=0;n<NfaceVertices;++n){
          int nid = faceVertices[f*NfaceVertices+n];
          sendFaces[id].v[n] = EToV[e*Nverts + nid];
        }

        std::sort(sendFaces[id].v, sendFaces[id].v+NfaceVertices,
                  std::less<hlong>());

        sendFaces[id].rank = rank;

        sendFaces[id].elementN = -1;
        sendFaces[id].faceN = -1;
        sendFaces[id].rankN = -1;

        ++Nsend[destRank];
      }
    }
  }

  // exchange byte counts
  comm.Alltoall(Nsend, Nrecv);

  // count incoming faces
  int allNrecv = 0;
  for(int rr=0;rr<size;++rr)
    allNrecv += Nrecv[rr];

  // find offsets for recv data
  for(int rr=1;rr<size;++rr)
    recvOffsets[rr] = recvOffsets[rr-1] + Nrecv[rr-1]; // byte offsets

  // buffer for incoming face data
  memory<face_t> recvFaces(allNrecv);

  // exchange parallel faces
  comm.Alltoallv(sendFaces, Nsend, sendOffsets,
                 recvFaces, Nrecv, recvOffsets);

  // local sort allNrecv received faces
  sort(recvFaces.ptr(), recvFaces.ptr()+allNrecv,
      [&](const face_t& a, const face_t& b) {
        return std::lexicographical_compare(a.v, a.v+NfaceVertices,
                                            b.v, b.v+NfaceVertices);
      });

  // find matches
  #pragma omp parallel for
  for(int n=0;n<allNrecv-1;++n){
    // since vertices are ordered we just look for pairs
    if(std::equal(recvFaces[n].v, recvFaces[n].v+NfaceVertices,
                  recvFaces[n+1].v)){
      recvFaces[n].elementN = recvFaces[n+1].element;
      recvFaces[n].faceN = recvFaces[n+1].face;
      recvFaces[n].rankN = recvFaces[n+1].rank;

      recvFaces[n+1].elementN = recvFaces[n].element;
      recvFaces[n+1].faceN = recvFaces[n].face;
      recvFaces[n+1].rankN = recvFaces[n].rank;
    }
  }

  // sort back to original ordering
  sort(recvFaces.ptr(), recvFaces.ptr()+allNrecv,
      [](const face_t& a, const face_t& b) {
        if(a.rank < b.rank) return true;
        if(a.rank > b.rank) return false;

        if(a.element < b.element) return true;
        if(a.element > b.element) return false;

        return (a.face < b.face);
      });

  // send faces back from whence they came
  comm.Alltoallv(recvFaces, Nrecv, recvOffsets,
                 sendFaces, Nsend, sendOffsets);

  // extract connectivity info
  #pragma omp parallel for
  for(dlong n=0;n<Nelements*Nfaces;++n)
    EToP[n] = -1;

  #pragma omp parallel for
  for(int n=0;n<allNsend;++n){
    dlong e = sendFaces[n].element;
    dlong eN = sendFaces[n].elementN;
    int f = sendFaces[n].face;
    int fN = sendFaces[n].faceN;
    int rN = sendFaces[n].rankN;

    if(e>=0 && f>=0 && eN>=0 && fN>=0){
      EToE[e*Nfaces+f] = eN;
      EToF[e*Nfaces+f] = fN;
      EToP[e*Nfaces+f] = rN;
    }
  }

  //record the number of elements in the whole mesh
  NelementsGlobal = Nelements;
  comm.Allreduce(NelementsGlobal);
}

} //namespace libp
