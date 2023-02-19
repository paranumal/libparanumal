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
void mesh_t::SurfaceGeometricFactorsHex3D(){

  /* unified storage array for geometric factors */
  Nsgeo = 8;

  NXID  = 0;
  NYID  = 1;
  NZID  = 2;
  SJID  = 3;
  IJID  = 4;
  IHID  = 5;
  WSJID = 6;
  WIJID = 7;

  props["defines/" "p_Nsgeo"]= Nsgeo;
  props["defines/" "p_NXID"]= NXID;
  props["defines/" "p_NYID"]= NYID;
  props["defines/" "p_NZID"]= NZID;
  props["defines/" "p_SJID"]= SJID;
  props["defines/" "p_IJID"]= IJID;
  props["defines/" "p_IHID"]= IHID;
  props["defines/" "p_WSJID"]= WSJID;
  props["defines/" "p_WIJID"]= WIJID;

  sgeo.malloc(Nelements*Nsgeo*Nfp*Nfaces);

  memory<dfloat> h((Nelements+totalHaloPairs)*Nfp*Nfaces);

  memory<dfloat> xre(Np);
  memory<dfloat> xse(Np);
  memory<dfloat> xte(Np);
  memory<dfloat> yre(Np);
  memory<dfloat> yse(Np);
  memory<dfloat> yte(Np);
  memory<dfloat> zre(Np);
  memory<dfloat> zse(Np);
  memory<dfloat> zte(Np);

  for(dlong e=0;e<Nelements;++e){ /* for each element */

    for(int k=0;k<Nq;++k){
      for(int j=0;j<Nq;++j){
        for(int i=0;i<Nq;++i){

          int n = i + j*Nq + k*Nq*Nq;
          xre[n] = 0; xse[n] = 0; xte[n] = 0;
          yre[n] = 0; yse[n] = 0; yte[n] = 0;
          zre[n] = 0; zse[n] = 0; zte[n] = 0;

          for(int m=0;m<Nq;++m){
            int idr = e*Np + k*Nq*Nq + j*Nq + m;
            int ids = e*Np + k*Nq*Nq + m*Nq + i;
            int idt = e*Np + m*Nq*Nq + j*Nq + i;
            xre[n] += D[i*Nq+m]*x[idr];
            xse[n] += D[j*Nq+m]*x[ids];
            xte[n] += D[k*Nq+m]*x[idt];
            yre[n] += D[i*Nq+m]*y[idr];
            yse[n] += D[j*Nq+m]*y[ids];
            yte[n] += D[k*Nq+m]*y[idt];
            zre[n] += D[i*Nq+m]*z[idr];
            zse[n] += D[j*Nq+m]*z[ids];
            zte[n] += D[k*Nq+m]*z[idt];
          }
        }
      }
    }

    for(int f=0;f<Nfaces;++f){ // for each face

      for(int i=0;i<Nfp;++i){  // for each node on face

        /* volume index of face node */
        int n = faceNodes[f*Nfp+i];

        dfloat xr = xre[n], xs = xse[n], xt = xte[n];
        dfloat yr = yre[n], ys = yse[n], yt = yte[n];
        dfloat zr = zre[n], zs = zse[n], zt = zte[n];

        /* determinant of Jacobian matrix */
        dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);

        dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
        dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
        dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;

        /* face f normal and length */
        dfloat nx=0.0, ny=0.0, nz=0.0;
        switch(f){
        case 0: nx = -tx; ny = -ty; nz = -tz; break;
        case 1: nx = -sx; ny = -sy; nz = -sz; break;
        case 2: nx = +rx; ny = +ry; nz = +rz; break;
        case 3: nx = +sx; ny = +sy; nz = +sz; break;
        case 4: nx = -rx; ny = -ry; nz = -rz; break;
        case 5: nx = +tx; ny = +ty; nz = +tz; break;
        }

        dfloat sJ = sqrt(nx*nx+ny*ny+nz*nz);
        nx /= sJ; ny /= sJ; nz /= sJ;
        sJ *= J;

        /* output index */
        dlong base = Nsgeo*(Nfaces*Nfp*e + Nfp*f + i);

        /* store normal, surface Jacobian, and reciprocal of volume Jacobian */
        sgeo[base+NXID] = nx;
        sgeo[base+NYID] = ny;
        sgeo[base+NZID] = nz;
        sgeo[base+SJID] = sJ;
        sgeo[base+IJID] = 1./J;

        sgeo[base+WIJID] = 1./(J*gllw[0]);
        sgeo[base+WSJID] = sJ*gllw[i%Nq]*gllw[i/Nq];

        h[Nfaces*Nfp*e + Nfp*f + i] = sJ/J;

        // computeFrame(nx, ny, nz,
        //              sgeo[base+STXID], sgeo[base+STYID], sgeo[base+STZID],
        //              sgeo[base+SBXID], sgeo[base+SBYID], sgeo[base+SBZID]);
      }
    }
  }

  halo.Exchange(h, Nfp*Nfaces);

  for(dlong e=0;e<Nelements;++e){ /* for each non-halo element */
    for(int f=0;f<Nfaces;++f){
      for(int n=0;n<Nfp;++n){
        dlong baseM = e*Nfp*Nfaces + f*Nfp + n;
        dlong baseP = mapP[baseM];
        if(baseP<0) baseP = baseM;

        // rescaling - missing factor of 2 ? (only impacts penalty and thus stiffness)
        dfloat hinvM = h[baseM];
        dfloat hinvP = h[baseP];
        sgeo[baseM*Nsgeo+IHID] = std::max(hinvM,hinvP);

        // if (EToB[f+e*Nfaces] > 0) { //enforce a stronger penalty on boundaries
        //   sgeo[baseM*Nsgeo+IHID] *= 2;
        // }
      }
    }
  }

  o_sgeo = platform.malloc<dfloat>(sgeo);

  if constexpr (std::is_same_v<dfloat,pfloat>) {
    o_pfloat_sgeo = o_sgeo;
  } else {
    memory<pfloat> pfloat_sgeo(Nelements*Nsgeo*Nfaces*Nfp);

    for(int n=0;n<Nsgeo*Nelements*Nfaces*Nfp;++n)
      pfloat_sgeo[n] = sgeo[n];

    o_pfloat_sgeo = platform.malloc<pfloat>(pfloat_sgeo);
  }
}

} //namespace libp
