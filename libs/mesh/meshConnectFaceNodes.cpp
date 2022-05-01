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

namespace libp {

// serial face-node to face-node connection
void mesh_t::ConnectFaceNodes(){

  /* Build the permutation array R */
  memory<int> R;

  switch (elementType) {
    case Mesh::TRIANGLES:
      FaceNodeMatchingTri2D(r, s, faceNodes, faceVertices, R);
      break;
    case Mesh::QUADRILATERALS:
      FaceNodeMatchingQuad2D(r, s, faceNodes, faceVertices, R);
      break;
    case Mesh::TETRAHEDRA:
      FaceNodeMatchingTet3D(r, s, t, faceNodes, faceVertices, R);
      break;
    case Mesh::HEXAHEDRA:
      FaceNodeMatchingHex3D(r, s, t, faceNodes, faceVertices, R);
      break;
  }

  /* volume indices of the interior and exterior face nodes for each element */
  vmapM.malloc(Nfp*Nfaces*Nelements);
  vmapP.malloc(Nfp*Nfaces*Nelements);
  mapP.malloc(Nfp*Nfaces*Nelements);

  /* assume elements already connected */
  #pragma omp parallel for collapse(2)
  for(dlong eM=0;eM<Nelements;++eM){
    for(int fM=0;fM<Nfaces;++fM){
      dlong eP = EToE[eM*Nfaces+fM];
      int fP = EToF[eM*Nfaces+fM];
      if(eP<0 || fP<0){ // fake connections for unconnected faces
        for(int nM=0;nM<Nfp;++nM){
          const int idM = faceNodes[fM*Nfp+nM];
          const dlong id = eM*Nfaces*Nfp + fM*Nfp + nM;
          vmapM[id] = idM + eM*Np;
          vmapP[id] = idM + eM*Np;
          mapP[id]  = id;
        }
      } else {
        //Find the rotation of the face from where the first vertex of the face is
        hlong vf0P = EToV[eP*Nverts + faceVertices[fP*NfaceVertices+0]];
        int rot=0;
        for (;rot<NfaceVertices;++rot) {
          hlong vfM = EToV[eM*Nverts + faceVertices[fM*NfaceVertices+rot]];
          if (vfM == vf0P) break;
        }

        /* for each node on this face use the permuation array
           to select the neighbor node */
        for(int nM=0;nM<Nfp;++nM){
          const int nP  = R[fM*Nfaces*Nfp*NfaceVertices
                            + fP*Nfp*NfaceVertices
                            + rot*Nfp + nM];
          const int idM = faceNodes[fM*Nfp+nM];
          const int idP = faceNodes[fP*Nfp+nP];

          const dlong id = eM*Nfaces*Nfp + fM*Nfp + nM;
          vmapM[id] = idM + eM*Np;
          vmapP[id] = idP + eP*Np;
          mapP[id]  = eP*Nfaces*Nfp + fP*Nfp + nP;
        }
      }
    }
  }

  o_vmapM = platform.malloc<dlong>(vmapM);
  o_vmapP = platform.malloc<dlong>(vmapP);
  o_mapP  = platform.malloc<dlong>(mapP);
}

} //namespace libp
