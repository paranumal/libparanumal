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

/* compute outwards facing normals, surface Jacobian, and volume Jacobian for all face nodes */
void mesh_t::SurfaceGeometricFactorsQuad2D(){

  /* unified storage array for geometric factors */
  Nsgeo = 7;

  NXID  = 0;
  NYID  = 1;
  SJID  = 2;
  IJID  = 3;
  IHID  = 4;
  WSJID = 5;
  WIJID = 6;

  props["defines/" "p_Nsgeo"]= Nsgeo;
  props["defines/" "p_NXID"]= NXID;
  props["defines/" "p_NYID"]= NYID;
  props["defines/" "p_SJID"]= SJID;
  props["defines/" "p_IJID"]= IJID;
  props["defines/" "p_IHID"]= IHID;
  props["defines/" "p_WSJID"]= WSJID;
  props["defines/" "p_WIJID"]= WIJID;

  sgeo.malloc(Nelements*Nsgeo*Nfp*Nfaces);

  memory<dfloat> hinv((Nelements+totalHaloPairs)*Nfp*Nfaces);

  memory<dfloat> xre(Np);
  memory<dfloat> xse(Np);
  memory<dfloat> yre(Np);
  memory<dfloat> yse(Np);

  for(dlong e=0;e<Nelements;++e){ /* for each element */

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

        hinv[Nfaces*Nfp*e + Nfp*f + i] = sJ/J;
      }
    }
  }

  halo.Exchange(hinv, Nfp*Nfaces);

  for(dlong e=0;e<Nelements;++e){ /* for each non-halo element */
    for(int f=0;f<Nfaces;++f){
      for(int n=0;n<Nfp;++n){
        dlong baseM = e*Nfp*Nfaces + f*Nfp + n;
        dlong baseP = mapP[baseM];
        if(baseP<0) baseP = baseM;

        // rescaling - missing factor of 2 ? (only impacts penalty and thus stiffness)
        dfloat hinvM = hinv[baseM];
        dfloat hinvP = hinv[baseP];
        sgeo[baseM*Nsgeo+IHID] = std::max(hinvM,hinvP);

        // if (EToB[f+e*Nfaces] > 0) { //enforce a stronger penalty on boundaries
        //   sgeo[baseM*Nsgeo+IHID] *= 2;
        // }
      }
    }
  }

  o_sgeo = platform.malloc<dfloat>(sgeo);
}

} //namespace libp
