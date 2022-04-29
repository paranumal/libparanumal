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

// custom geometric factors specialized for 3D tri on sphere

namespace libp {

void mesh_t::GeometricFactorsTri3D(){

  /*Set offsets*/
  Nvgeo = 10;

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

  #pragma omp parallel for
  for(int e=0;e<Nelements;++e){ /* for each element */

    for(int n=0;n<Np;++n){

      dfloat xn = x[n+e*Np];
      dfloat yn = y[n+e*Np];
      dfloat zn = z[n+e*Np];

      dfloat xr = 0, yr = 0, zr = 0;
      dfloat xs = 0, ys = 0, zs = 0;

      for(int m=0;m<Np;++m){

        dfloat Drnm = Dr[n*Np+m];
        dfloat Dsnm = Ds[n*Np+m];

        xr += Drnm*x[m+e*Np];
        yr += Drnm*y[m+e*Np];
        zr += Drnm*z[m+e*Np];

        xs += Dsnm*x[m+e*Np];
        ys += Dsnm*y[m+e*Np];
        zs += Dsnm*z[m+e*Np];

      }

      dfloat rx = ys*zn - zs*yn; // dXds x X
      dfloat ry = zs*xn - xs*zn;
      dfloat rz = xs*yn - ys*xn;

      dfloat sx = zr*yn - yr*zn; // -dXdr x X
      dfloat sy = xr*zn - zr*xn;
      dfloat sz = yr*xn - xr*yn;

      dfloat tx = yr*zs - zr*ys; // dXdr x dXds ~ X*|dXdr x dXds|/|X|
      dfloat ty = zr*xs - xr*zs;
      dfloat tz = xr*ys - yr*xs;

      dfloat Gx = tx, Gy = ty, Gz = tz;

      dfloat J = xn*tx + yn*ty + zn*tz;

      LIBP_ABORT("Negative J found at element " << e, J<1e-8);

      rx /= J;
      ry /= J;
      rz /= J;

      sx /= J;
      sy /= J;
      sz /= J;

      tx /= J;
      ty /= J;
      tz /= J;

      // use this for "volume" Jacobian
      J = sqrt(Gx*Gx+Gy*Gy+Gz*Gz);

      LIBP_ABORT("Negative J found at element " << e, J<1e-8);

      /* store geometric factors */
      int base = Nvgeo*Np*e + n;

      vgeo[base + Np*RXID] = rx;
      vgeo[base + Np*RYID] = ry;
      vgeo[base + Np*RZID] = rz;
      vgeo[base + Np*SXID] = sx;
      vgeo[base + Np*SYID] = sy;
      vgeo[base + Np*SZID] = sz;
      vgeo[base + Np*TXID] = tx;
      vgeo[base + Np*TYID] = ty;
      vgeo[base + Np*TZID] = tz;
      vgeo[base + Np*JID]  = J;

      ggeo[Nggeo*e*Np+n + Np*G00ID] = J*(rx*rx + ry*ry + rz*rz);
      ggeo[Nggeo*e*Np+n + Np*G01ID] = J*(rx*sx + ry*sy + rz*sz);
      ggeo[Nggeo*e*Np+n + Np*G11ID] = J*(sx*sx + sy*sy + sz*sz);

      wJ[e*Np + n]  = J;
    }
  }

  halo.Exchange(vgeo, Nvgeo*Np);

  o_wJ   = platform.malloc<dfloat>(wJ);
  o_vgeo = platform.malloc<dfloat>(vgeo);
  o_ggeo = platform.malloc<dfloat>(ggeo);
}

} //namespace libp
