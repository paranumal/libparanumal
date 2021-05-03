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

// serial face-vertex to face-vertex connection
void mesh_t::ConnectFaceVertices(){

  //allocate and fill a halo region in element-to-vertex mapping
  EToV = (hlong*) realloc(EToV, (Nelements+totalHaloPairs)*Nverts*sizeof(hlong));
  halo->Exchange(EToV, Nverts, ogs_hlong);

  /* volume indices of the interior and exterior face vertices for each element */
  VmapM = (dlong*) malloc(NfaceVertices*Nfaces*Nelements*sizeof(dlong));
  VmapP = (dlong*) malloc(NfaceVertices*Nfaces*Nelements*sizeof(dlong));

  /* assume elements already connected */
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      dlong eP = EToE[e*Nfaces+f];
      int fP = EToF[e*Nfaces+f];
      if(eP<0 || fP<0){ // fake connections for unconnected faces
        eP = e;
        fP = f;
      }

      /* for each vertex on this face find the neighbor vertex */
      for(int n=0;n<NfaceVertices;++n){
        dlong idM = faceVertices[f*NfaceVertices+n] + e*Nverts;
        hlong vM  = EToV[idM];

        dlong idP=idM;
        for(int m=0;m<NfaceVertices;++m){
          idP = faceVertices[fP*NfaceVertices+m] + eP*Nverts;
          if (EToV[idP]==vM) break;
        }

        dlong id = Nfaces*NfaceVertices*e + f*NfaceVertices + n;
        VmapM[id] = idM;
        VmapP[id] = idP;
      }
    }
  }
}

