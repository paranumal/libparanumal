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

static int findBestMatch(dfloat x1, dfloat y1, dfloat z1,
                   int Np2, int *nodeList, dfloat *x2, dfloat *y2, dfloat *z2, int *nP){

  int matchIndex;
  dfloat mindist2=1e9;

  for(int n=0;n<Np2;++n){

    /* next node */
    const int i2 = nodeList[n];

    /* distance between target and next node */
    const dfloat dist2 = pow(x1-x2[i2],2) + pow(y1-y2[i2],2) + pow(z1-z2[i2],2);

    /* if next node is closer to target update match */
    if(n==0 || dist2<mindist2){
      mindist2 = dist2;
      matchIndex = i2;
      *nP = n;
    }
  }
  if(mindist2>1e-3) {
    stringstream ss;
    ss << "Bad match: x,y,z = " << x1 << ", " << y1 << ", " << z1 << "\n";
    LIBP_ABORT(ss.str())
  }

  return matchIndex;
}


// serial face-node to face-node connection
void mesh3D::ConnectFaceNodes(){

  /* volume indices of the interior and exterior face nodes for each element */
  vmapM = (dlong*) calloc(Nfp*Nfaces*Nelements, sizeof(dlong));
  vmapP = (dlong*) calloc(Nfp*Nfaces*Nelements, sizeof(dlong));
  mapP  = (dlong*) calloc(Nfp*Nfaces*Nelements, sizeof(dlong));

  /* assume elements already connected */
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      dlong eP = EToE[e*Nfaces+f];
      int fP = EToF[e*Nfaces+f];
      if(eP<0 || fP<0){ // fake connections for unconnected faces
        eP = e;
        fP = f;
      }
      /* for each node on this face find the neighbor node */
      for(int n=0;n<Nfp;++n){
        dlong  idM = faceNodes[f*Nfp+n] + e*Np;
        dfloat xM = x[idM];
        dfloat yM = y[idM];
        dfloat zM = z[idM];
        int nP;

        int  idP = findBestMatch(xM, yM, zM,
				 Nfp,
				 faceNodes+fP*Nfp,
				 x+eP*Np,
				 y+eP*Np,
				 z+eP*Np, &nP);

        dlong id = Nfaces*Nfp*e + f*Nfp + n;
        vmapM[id] = idM;
        vmapP[id] = idP + eP*Np;
        mapP[id] = eP*Nfaces*Nfp + fP*Nfp + nP;
      }
    }
  }
}

