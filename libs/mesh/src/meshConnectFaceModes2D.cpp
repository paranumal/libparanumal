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
#include "mesh2D.hpp"

// serial face-mode to face-mode connection
void mesh2D::ConnectFaceModes(int *faceModes, dfloat *V){

  /* volume indices of the interior and exterior face modes for each element */
  mmapM = (dlong*) calloc(Nfp*Nfaces*Nelements, sizeof(dlong));
  mmapP = (dlong*) calloc(Nfp*Nfaces*Nelements, sizeof(dlong));
  mmapS = (int*) calloc(Nfp*Nfaces*Nelements, sizeof(int));

  dfloat *VM = (dfloat *) calloc(Np,sizeof(dfloat));
  dfloat *VP = (dfloat *) calloc(Np,sizeof(dfloat));

  /* assume elements already connected */
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      dlong eP = EToE[e*Nfaces+f];
      int fP = EToF[e*Nfaces+f];
      if(eP<0 || fP<0){ // fake connections for unconnected faces
        eP = e;
        fP = f;
      }

      /* for each mode on this face find the neighbor mode */
      for(int n=0;n<Nfp;++n){
        int m = faceModes[n+f*Nfp]; //get face mode number

        for (int i=0;i<Nfp;i++) {
          int k = faceNodes[i+f*Nfp];
          VM[i] = V[m+k*Np]; //evaluate mode at WB nodes on face
        }

        dfloat mindist = 1E9;
        int sflag=0;
        int mMatch=0;
        for (int nP=0;nP<Nfp;nP++) {
          //test the modes on face fP
          int mP = faceModes[nP+fP*Nfp];

          for (int i=0;i<Nfp;i++) {
            //get neighbouring node
            dlong id = i+f*Nfp+e*Nfp*Nfaces;
            int k = vmapP[id]%Np;

            VP[i] = V[mP+k*Np]; //evaluate mode at WB nodes on face
          }

          dfloat dist1=0, dist2=0;
          for (int i=0;i<Nfp;i++){
            dist1 += pow(VM[i]-VP[i],2);
            dist2 += pow(VM[i]+VP[i],2);
          }
          dist1 = sqrt(dist1);
          dist2 = sqrt(dist2);

          /* if next node is closer to target update match */
          if(dist1<mindist){
            mindist = dist1;
            mMatch = mP;
            sflag = 1;
          }
          if(dist2<mindist){
            mindist = dist2;
            mMatch = mP;
            sflag = -1;
          }
        }
        if(mindist>1e-3) {
          stringstream ss;
          ss << "Bad match: e = " << e
             << ", f = "<< f
             << ", mode = " << m << "\n";
          LIBP_ABORT(ss.str())
        }

        dlong id  = Nfaces*Nfp*e + f*Nfp + n;
        dlong idM = faceModes[f*Nfp+n] + e*Np;
        dlong idP = mMatch + eP*Np;

        mmapM[id] = idM;
        mmapP[id] = idP;
        mmapS[id] = sflag;
      }
    }
  }
}
