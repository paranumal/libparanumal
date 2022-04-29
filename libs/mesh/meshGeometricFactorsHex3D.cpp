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

void mesh_t::GeometricFactorsHex3D(){

  /*Set offsets*/
  Nvgeo = 12;

  RXID  = 0;
  RYID  = 1;
  RZID  = 2;
  SXID  = 3;
  SYID  = 4;
  SZID  = 5;
  TXID  = 6;
  TYID  = 7;
  TZID  = 8;
  JID   = 9;
  JWID  = 10;
  IJWID = 11;

  props["defines/" "p_Nvgeo"]= Nvgeo;
  props["defines/" "p_RXID"]= RXID;
  props["defines/" "p_SXID"]= SXID;
  props["defines/" "p_TXID"]= TXID;

  props["defines/" "p_RYID"]= RYID;
  props["defines/" "p_SYID"]= SYID;
  props["defines/" "p_TYID"]= TYID;

  props["defines/" "p_RZID"]= RZID;
  props["defines/" "p_SZID"]= SZID;
  props["defines/" "p_TZID"]= TZID;

  props["defines/" "p_JID"]= JID;
  props["defines/" "p_JWID"]= JWID;
  props["defines/" "p_IJWID"]= IJWID;

  /* unified storage array for geometric factors */
  /* note that we have volume geometric factors for each node */
  vgeo.malloc((Nelements+totalHaloPairs)*Nvgeo*Np);

  Nggeo = 6;

  G00ID=0;
  G01ID=1;
  G02ID=2;
  G11ID=3;
  G12ID=4;
  G22ID=5;

  props["defines/" "p_Nggeo"]= Nggeo;
  props["defines/" "p_G00ID"]= G00ID;
  props["defines/" "p_G01ID"]= G01ID;
  props["defines/" "p_G02ID"]= G02ID;
  props["defines/" "p_G11ID"]= G11ID;
  props["defines/" "p_G12ID"]= G12ID;
  props["defines/" "p_G22ID"]= G22ID;

  /* number of second order geometric factors */
  ggeo.malloc(Nelements*Nggeo*Np);

  wJ.malloc(Nelements*Np);

  // dfloat minJ = 1e9, maxJ = -1e9, maxSkew = 0;

  #pragma omp parallel for
  for(dlong e=0;e<Nelements;++e){ /* for each element */

    for(int k=0;k<Nq;++k){
      for(int j=0;j<Nq;++j){
        for(int i=0;i<Nq;++i){

          int n = i + j*Nq + k*Nq*Nq;

          dfloat xr = 0, xs = 0, xt = 0;
          dfloat yr = 0, ys = 0, yt = 0;
          dfloat zr = 0, zs = 0, zt = 0;
          for(int m=0;m<Nq;++m){
            int idr = e*Np + k*Nq*Nq + j*Nq + m;
            int ids = e*Np + k*Nq*Nq + m*Nq + i;
            int idt = e*Np + m*Nq*Nq + j*Nq + i;
            xr += D[i*Nq+m]*x[idr];
            xs += D[j*Nq+m]*x[ids];
            xt += D[k*Nq+m]*x[idt];
            yr += D[i*Nq+m]*y[idr];
            ys += D[j*Nq+m]*y[ids];
            yt += D[k*Nq+m]*y[idt];
            zr += D[i*Nq+m]*z[idr];
            zs += D[j*Nq+m]*z[ids];
            zt += D[k*Nq+m]*z[idt];
          }

          /* compute geometric factors for affine coordinate transform*/
          dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);

          // dfloat hr = sqrt(xr*xr+yr*yr+zr*zr);
          // dfloat hs = sqrt(xs*xs+ys*ys+zs*zs);
          // dfloat ht = sqrt(xt*xt+yt*yt+zt*zt);
          // minJ = mymin(J, minJ);
          // maxJ = mymax(J, maxJ);
          // maxSkew = mymax(maxSkew, hr/hs);
          // maxSkew = mymax(maxSkew, hr/ht);
          // maxSkew = mymax(maxSkew, hs/hr);
          // maxSkew = mymax(maxSkew, hs/ht);
          // maxSkew = mymax(maxSkew, ht/hr);
          // maxSkew = mymax(maxSkew, ht/hs);

          LIBP_ABORT("Negative J found at element " << e,
                     J<1e-12);

          dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
          dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
          dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;

          dfloat JW = J*gllw[i]*gllw[j]*gllw[k];

          /* store geometric factors */
          vgeo[Nvgeo*Np*e + n + Np*RXID] = rx;
          vgeo[Nvgeo*Np*e + n + Np*RYID] = ry;
          vgeo[Nvgeo*Np*e + n + Np*RZID] = rz;

          vgeo[Nvgeo*Np*e + n + Np*SXID] = sx;
          vgeo[Nvgeo*Np*e + n + Np*SYID] = sy;
          vgeo[Nvgeo*Np*e + n + Np*SZID] = sz;

          vgeo[Nvgeo*Np*e + n + Np*TXID] = tx;
          vgeo[Nvgeo*Np*e + n + Np*TYID] = ty;
          vgeo[Nvgeo*Np*e + n + Np*TZID] = tz;

          vgeo[Nvgeo*Np*e + n + Np*JID]  = J;
          vgeo[Nvgeo*Np*e + n + Np*JWID] = JW;
          vgeo[Nvgeo*Np*e + n + Np*IJWID] = 1./JW;

          /* store second order geometric factors */
          ggeo[Nggeo*Np*e + n + Np*G00ID] = JW*(rx*rx + ry*ry + rz*rz);
          ggeo[Nggeo*Np*e + n + Np*G01ID] = JW*(rx*sx + ry*sy + rz*sz);
          ggeo[Nggeo*Np*e + n + Np*G02ID] = JW*(rx*tx + ry*ty + rz*tz);
          ggeo[Nggeo*Np*e + n + Np*G11ID] = JW*(sx*sx + sy*sy + sz*sz);
          ggeo[Nggeo*Np*e + n + Np*G12ID] = JW*(sx*tx + sy*ty + sz*tz);
          ggeo[Nggeo*Np*e + n + Np*G22ID] = JW*(tx*tx + ty*ty + tz*tz);

          wJ[Np*e + n] = JW;
        }
      }
    }
  }

  halo.Exchange(vgeo, Nvgeo*Np);

  o_wJ   = platform.malloc<dfloat>(wJ);
  o_vgeo = platform.malloc<dfloat>(vgeo);
  o_ggeo = platform.malloc<dfloat>(ggeo);


  #if 0
    dfloat globalMinJ, globalMaxJ, globalMaxSkew;

    MPI_Reduce(&minJ, &globalMinJ, 1, MPI_DFLOAT, MPI_MIN, 0, comm);
    MPI_Reduce(&maxJ, &globalMaxJ, 1, MPI_DFLOAT, MPI_MAX, 0, comm);
    MPI_Reduce(&maxSkew, &globalMaxSkew, 1, MPI_DFLOAT, MPI_MAX, 0, comm);

    if(rank==0)
      printf("J in range [%g,%g] and max Skew = %g\n", globalMinJ, globalMaxJ, globalMaxSkew);
  #endif
}

} //namespace libp
