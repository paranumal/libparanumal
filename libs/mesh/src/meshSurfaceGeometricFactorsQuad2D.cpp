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

/* compute outwards facing normals, surface Jacobian, and volume Jacobian for all face nodes */
void meshQuad2D::SurfaceGeometricFactors(){

  /* unified storage array for geometric factors */
  Nsgeo = 7;
  sgeo = (dfloat*) calloc((Nelements+totalHaloPairs)*
                                Nsgeo*Nfp*Nfaces,
                                sizeof(dfloat));

  dfloat *xre = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *xse = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *yre = (dfloat*) calloc(Np, sizeof(dfloat));
  dfloat *yse = (dfloat*) calloc(Np, sizeof(dfloat));

  for(dlong e=0;e<Nelements+totalHaloPairs;++e){ /* for each element */

    for(int j=0;j<Nq;++j){
      for(int i=0;i<Nq;++i){
        int n = i + j*Nq;
        xre[n] = 0; xse[n] = 0;
        yre[n] = 0; yse[n] = 0;

        for(int m=0;m<Nq;++m){
          int idr = e*Np + j*Nq + m;
          int ids = e*Np + m*Nq + i;
          xre[n] += D[i*Nq+m]*x[idr];
          xse[n] += D[j*Nq+m]*x[ids];
          yre[n] += D[i*Nq+m]*y[idr];
          yse[n] += D[j*Nq+m]*y[ids];
        }
      }
    }

    for(int f=0;f<Nfaces;++f){ // for each face
      for(int i=0;i<Nq;++i){  // for each node on face

        /* volume index of face node */
        int n = faceNodes[f*Nfp+i];

        dfloat xr = xre[n], xs = xse[n];
        dfloat yr = yre[n], ys = yse[n];

        /* compute geometric factors for affine coordinate transform*/
        dfloat J = xr*ys - xs*yr;

        dfloat rx =  ys/J;
        dfloat ry = -xs/J;
        dfloat sx = -yr/J;
        dfloat sy =  xr/J;

        /* face f normal and length */
        dfloat nx=0.0, ny=0.0;
        switch(f){
        case 0: nx = -sx; ny = -sy; break;
        case 1: nx = +rx; ny = +ry; break;
        case 2: nx = +sx; ny = +sy; break;
        case 3: nx = -rx; ny = -ry; break;
        }
        dfloat  sJ = sqrt((nx)*(nx)+(ny)*(ny));
        nx /= sJ; ny /= sJ;
        sJ *= J;

        /* output index */
        dlong base = Nsgeo*(Nfaces*Nfp*e + Nfp*f + i);

        /* store normal, surface Jacobian, and reciprocal of volume Jacobian */
        sgeo[base+NXID] = nx;
        sgeo[base+NYID] = ny;
        sgeo[base+SJID] = sJ;
        sgeo[base+IJID] = 1./J;

        sgeo[base+WIJID] = 1./(J*gllw[0]);
        sgeo[base+WSJID] = sJ*gllw[i];
      }
    }
  }

  for(dlong e=0;e<Nelements;++e){ /* for each non-halo element */
    for(int n=0;n<Nfp*Nfaces;++n){
      dlong baseM = e*Nfp*Nfaces + n;
      dlong baseP = mapP[baseM];
      if(baseP<0) baseP = baseM;

      // rescaling - missing factor of 2 ? (only impacts penalty and thus stiffness)
      dfloat hinvM = sgeo[baseM*Nsgeo + SJID]*sgeo[baseM*Nsgeo + IJID];
      dfloat hinvP = sgeo[baseP*Nsgeo + SJID]*sgeo[baseP*Nsgeo + IJID];

      sgeo[baseM*Nsgeo+IHID] = mymax(hinvM,hinvP);
      sgeo[baseP*Nsgeo+IHID] = mymax(hinvM,hinvP);
    }
  }
}
