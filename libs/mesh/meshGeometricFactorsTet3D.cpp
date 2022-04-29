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

void mesh_t::GeometricFactorsTet3D(){

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
  vgeo.malloc((Nelements+totalHaloPairs)*Nvgeo);

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
  ggeo.malloc(Nelements*Nggeo);

  wJ.malloc(Nelements);


  // dfloat minJ = 1e9, maxJ = -1e9;

  #pragma omp parallel for
  for(dlong e=0;e<Nelements;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    dlong id = e*Nverts;

    /* vertex coordinates */
    dfloat xe1 = EX[id+0], ye1 = EY[id+0], ze1 = EZ[id+0];
    dfloat xe2 = EX[id+1], ye2 = EY[id+1], ze2 = EZ[id+1];
    dfloat xe3 = EX[id+2], ye3 = EY[id+2], ze3 = EZ[id+2];
    dfloat xe4 = EX[id+3], ye4 = EY[id+3], ze4 = EZ[id+3];

    /* Jacobian matrix */
    dfloat xr = 0.5*(xe2-xe1), xs = 0.5*(xe3-xe1), xt = 0.5*(xe4-xe1);
    dfloat yr = 0.5*(ye2-ye1), ys = 0.5*(ye3-ye1), yt = 0.5*(ye4-ye1);
    dfloat zr = 0.5*(ze2-ze1), zs = 0.5*(ze3-ze1), zt = 0.5*(ze4-ze1);

    /* compute geometric factors for affine coordinate transform*/
    dfloat J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt);

    dfloat rx =  (ys*zt - zs*yt)/J, ry = -(xs*zt - zs*xt)/J, rz =  (xs*yt - ys*xt)/J;
    dfloat sx = -(yr*zt - zr*yt)/J, sy =  (xr*zt - zr*xt)/J, sz = -(xr*yt - yr*xt)/J;
    dfloat tx =  (yr*zs - zr*ys)/J, ty = -(xr*zs - zr*xs)/J, tz =  (xr*ys - yr*xs)/J;

    LIBP_ABORT("Negative J found at element " << e, J<0);

    // minJ = mymin(minJ,J);
    // maxJ = mymax(maxJ,J);

    /* store geometric factors */
    vgeo[Nvgeo*e + RXID] = rx;
    vgeo[Nvgeo*e + RYID] = ry;
    vgeo[Nvgeo*e + RZID] = rz;
    vgeo[Nvgeo*e + SXID] = sx;
    vgeo[Nvgeo*e + SYID] = sy;
    vgeo[Nvgeo*e + SZID] = sz;
    vgeo[Nvgeo*e + TXID] = tx;
    vgeo[Nvgeo*e + TYID] = ty;
    vgeo[Nvgeo*e + TZID] = tz;
    vgeo[Nvgeo*e +  JID] = J;
    //    printf("geo: %g,%g,%g - %g,%g,%g - %g,%g,%g\n",
    //     rx,ry,rz, sx,sy,sz, tx,ty,tz);

    /* store second order geometric factors */
    ggeo[Nggeo*e + G00ID] = J*(rx*rx + ry*ry + rz*rz);
    ggeo[Nggeo*e + G01ID] = J*(rx*sx + ry*sy + rz*sz);
    ggeo[Nggeo*e + G02ID] = J*(rx*tx + ry*ty + rz*tz);
    ggeo[Nggeo*e + G11ID] = J*(sx*sx + sy*sy + sz*sz);
    ggeo[Nggeo*e + G12ID] = J*(sx*tx + sy*ty + sz*tz);
    ggeo[Nggeo*e + G22ID] = J*(tx*tx + ty*ty + tz*tz);

    wJ[e] = J;
  }
  //printf("minJ = %g, maxJ = %g\n", minJ, maxJ);

  halo.Exchange(vgeo, Nvgeo);

  o_wJ   = platform.malloc<dfloat>(wJ);
  o_vgeo = platform.malloc<dfloat>(vgeo);
  o_ggeo = platform.malloc<dfloat>(ggeo);
}

} //namespace libp
