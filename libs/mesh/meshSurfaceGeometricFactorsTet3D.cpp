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

void mesh_t::SurfaceGeometricFactorsTet3D(){

  /* unified storage array for geometric factors */
  Nsgeo = 6;

  NXID  = 0;
  NYID  = 1;
  NZID  = 2;
  SJID  = 3;
  IJID  = 4;
  IHID  = 5;

  props["defines/" "p_Nsgeo"]= Nsgeo;
  props["defines/" "p_NXID"]= NXID;
  props["defines/" "p_NYID"]= NYID;
  props["defines/" "p_NZID"]= NZID;
  props["defines/" "p_SJID"]= SJID;
  props["defines/" "p_IJID"]= IJID;
  props["defines/" "p_IHID"]= IHID;

  sgeo.malloc(Nelements*Nsgeo*Nfaces);

  memory<dfloat> hinv((Nelements+totalHaloPairs)*Nfaces);

  for(dlong e=0;e<Nelements;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    dlong id = e*Nverts;
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

    LIBP_ABORT("Negative J found at element " << e,
               J<0);

    /* face 1 */
    dlong base = Nsgeo*Nfaces*e;
    dfloat nx1 = -tx;
    dfloat ny1 = -ty;
    dfloat nz1 = -tz;
    dfloat sJ1 = sqrt((nx1)*(nx1)+(ny1)*(ny1)+(nz1)*(nz1));

    sgeo[base+NXID] = nx1/sJ1;
    sgeo[base+NYID] = ny1/sJ1;
    sgeo[base+NZID] = nz1/sJ1;
    sgeo[base+SJID] = sJ1*J;
    sgeo[base+IJID] = 1./J;

    hinv[Nfaces*e+0] = 0.5*sJ1;

    /* face 2 */
    base += Nsgeo;
    dfloat nx2 = -sx;
    dfloat ny2 = -sy;
    dfloat nz2 = -sz;
    dfloat sJ2 = sqrt((nx2)*(nx2)+(ny2)*(ny2)+(nz2)*(nz2));

    sgeo[base+NXID] = nx2/sJ2;
    sgeo[base+NYID] = ny2/sJ2;
    sgeo[base+NZID] = nz2/sJ2;
    sgeo[base+SJID] = sJ2*J;
    sgeo[base+IJID] = 1./J;

    hinv[Nfaces*e+1] = 0.5*sJ2;

    /* face 3 */
    base += Nsgeo;
    dfloat nx3 = rx+sx+tx;
    dfloat ny3 = ry+sy+ty;
    dfloat nz3 = rz+sz+tz;
    dfloat sJ3 = sqrt((nx3)*(nx3)+(ny3)*(ny3)+(nz3)*(nz3));

    sgeo[base+NXID] = nx3/sJ3;
    sgeo[base+NYID] = ny3/sJ3;
    sgeo[base+NZID] = nz3/sJ3;
    sgeo[base+SJID] = sJ3*J;
    sgeo[base+IJID] = 1./J;

    hinv[Nfaces*e+2] = 0.5*sJ3;

    /* face 4 */
    base += Nsgeo;
    dfloat nx4 = -rx;
    dfloat ny4 = -ry;
    dfloat nz4 = -rz;
    dfloat sJ4 = sqrt((nx4)*(nx4)+(ny4)*(ny4)+(nz4)*(nz4));

    sgeo[base+NXID] = nx4/sJ4;
    sgeo[base+NYID] = ny4/sJ4;
    sgeo[base+NZID] = nz4/sJ4;
    sgeo[base+SJID] = sJ4*J;
    sgeo[base+IJID] = 1./J;

    hinv[Nfaces*e+3] = 0.5*sJ4;

#if 0
    printf("N1=(%g,%g,%g),sJ1=%g\n", nx1/sJ1,ny1/sJ1,nz1/sJ1,sJ1*J);
    printf("N2=(%g,%g,%g),sJ2=%g\n", nx2/sJ2,ny2/sJ2,nz2/sJ2,sJ2*J);
    printf("N3=(%g,%g,%g),sJ3=%g\n", nx3/sJ3,ny3/sJ3,nz3/sJ3,sJ3*J);
    printf("N4=(%g,%g,%g),sJ4=%g\n", nx4/sJ4,ny4/sJ4,nz4/sJ4,sJ4*J);
#endif
  }

  halo.Exchange(hinv, Nfaces);

  for(dlong eM=0;eM<Nelements;++eM){ /* for each non-halo element */
    for(int fM=0;fM<Nfaces;++fM){
      dlong eP = EToE[eM*Nfaces+fM];

      if (eP<0) eP = eM;

      int fP = EToF[eM*Nfaces+fM];
      if (fP<0) fP = fM;

      dlong baseM = eM*Nfaces + fM;
      dlong baseP = eP*Nfaces + fP;

      // rescaling - A = L*h/2 => (J*2) = (sJ*2)*h/2 => h  = 2*J/sJ
      dfloat hinvM = hinv[baseM];
      dfloat hinvP = hinv[baseP];
      sgeo[baseM*Nsgeo+IHID] = std::max(hinvM,hinvP);

      // if (EToB[fM+eM*Nfaces] > 0) { //enforce a stronger penalty on boundaries
      //   sgeo[baseM*Nsgeo+IHID] *= 2;
      // }
    }
  }

  o_sgeo = platform.malloc<dfloat>(sgeo);
}

} //namespace libp
