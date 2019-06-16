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
#include "mesh3D.hpp"

// custom geometric factors specialized for 3D quad on sphere

void meshQuad3D::GeometricFactors(){

  /* unified storage array for geometric factors */
  Nvgeo = 12; //

  /* note that we have volume geometric factors for each node */
  vgeo = (dfloat*) calloc((Nelements+totalHaloPairs)*Nvgeo*Np, sizeof(dfloat));

  cubvgeo = (dfloat*) calloc(Nelements*Nvgeo*cubNp, sizeof(dfloat));

  // Can be computed on the fly
  Nggeo = 7;
  ggeo  = (dfloat *) calloc(Nelements*Np*Nggeo, sizeof(dfloat));

  dfloat *cxr = (dfloat*) calloc(cubNq*cubNq, sizeof(dfloat));
  dfloat *cxs = (dfloat*) calloc(cubNq*cubNq, sizeof(dfloat));
  dfloat *cyr = (dfloat*) calloc(cubNq*cubNq, sizeof(dfloat));
  dfloat *cys = (dfloat*) calloc(cubNq*cubNq, sizeof(dfloat));
  dfloat *czr = (dfloat*) calloc(cubNq*cubNq, sizeof(dfloat));
  dfloat *czs = (dfloat*) calloc(cubNq*cubNq, sizeof(dfloat));
  dfloat *cx  = (dfloat*) calloc(cubNq*cubNq, sizeof(dfloat));
  dfloat *cy  = (dfloat*) calloc(cubNq*cubNq, sizeof(dfloat));
  dfloat *cz  = (dfloat*) calloc(cubNq*cubNq, sizeof(dfloat));

  for(int e=0;e<Nelements;++e){ /* for each element */

    for(int n=0;n<cubNq*cubNq;++n){
      cxr[n] = 0; cyr[n] = 0; czr[n] = 0;
      cxs[n] = 0; cys[n] = 0; czs[n] = 0;
      cx[n] = 0;  cy[n] = 0;  cz[n] = 0;
    }

    for(int j=0;j<Nq;++j){
      for(int i=0;i<Nq;++i){

  dfloat xij = x[i+j*Nq+e*Np];
  dfloat yij = y[i+j*Nq+e*Np];
  dfloat zij = z[i+j*Nq+e*Np];

  dfloat xr = 0, yr = 0, zr = 0;
  dfloat xs = 0, ys = 0, zs = 0;

  for(int n=0;n<Nq;++n){

    dfloat Din = D[i*Nq+n];
    dfloat Djn = D[j*Nq+n];

    xr += Din*x[n+j*Nq+e*Np];
    yr += Din*y[n+j*Nq+e*Np];
    zr += Din*z[n+j*Nq+e*Np];

    xs += Djn*x[i+n*Nq+e*Np];
    ys += Djn*y[i+n*Nq+e*Np];
    zs += Djn*z[i+n*Nq+e*Np];

  }

  {
    dfloat rx = ys*zij - zs*yij; // dXds x X
    dfloat ry = zs*xij - xs*zij;
    dfloat rz = xs*yij - ys*xij;

    dfloat sx = zr*yij - yr*zij; // -dXdr x X
    dfloat sy = xr*zij - zr*xij;
    dfloat sz = yr*xij - xr*yij;

    dfloat tx = yr*zs - zr*ys; // dXdr x dXds ~ X*|dXdr x dXds|/|X|
    dfloat ty = zr*xs - xr*zs;
    dfloat tz = xr*ys - yr*xs;

    dfloat Gx = tx, Gy = ty, Gz = tz;

    dfloat J = xij*tx + yij*ty + zij*tz;

    if(J<1e-8) {
      stringstream ss;
      ss << "Negative J found at element " << e << "\n";
      LIBP_ABORT(ss.str())
    }

    rx /= J;      sx /= J;      tx /= J;
    ry /= J;      sy /= J;      ty /= J;
    rz /= J;      sz /= J;      tz /= J;

    // use this for "volume" Jacobian
    dfloat Jnew = sqrt(Gx*Gx+Gy*Gy+Gz*Gz);  //(difference between actual Jacobian and sphere Jac)
    J = Jnew;

    if(J<1e-8) {
      stringstream ss;
      ss << "Negative J found at element " << e << "\n";
      LIBP_ABORT(ss.str())
    }
    //    printf("before: grad r = %g,%g,%g\n", rx, ry, rz);
  }

  dfloat GG00 = xr*xr+yr*yr+zr*zr;
  dfloat GG11 = xs*xs+ys*ys+zs*zs;
  dfloat GG01 = xr*xs+yr*ys+zr*zs;
  dfloat detGG = GG00*GG11 - GG01*GG01;

  // are these tangential
  dfloat rx = (xr*GG11-xs*GG01)/detGG;
  dfloat ry = (yr*GG11-ys*GG01)/detGG;
  dfloat rz = (zr*GG11-zs*GG01)/detGG;

  dfloat sx = (-xr*GG01+xs*GG00)/detGG;
  dfloat sy = (-yr*GG01+ys*GG00)/detGG;
  dfloat sz = (-zr*GG01+zs*GG00)/detGG;

  dfloat tx = yr*zs - zr*ys; // dXdr x dXds ~ X*|dXdr x dXds|/|X|
  dfloat ty = zr*xs - xr*zs;
  dfloat tz = xr*ys - yr*xs;

  // use this for "volume" Jacobian
  dfloat J = sqrt(tx*tx+ty*ty+tz*tz); // (difference between actual Jacobian and sphere Jac)

  //  printf("after: grad r = %g,%g,%g\n", rx, ry, rz);

  dfloat JW = J*gllw[i]*gllw[j];

  /* store geometric factors */
  int base = Nvgeo*Np*e + j*Nq + i;

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
  vgeo[base + Np*JWID] = JW;
  vgeo[base + Np*IJWID] = 1./JW;

  /* store second order geometric factors (can be computed on the fly, later!!!)*/
  int gbase = Nggeo*Np*e + j*Nq + i;
  ggeo[gbase + Np*G00ID] = JW*(rx*rx + ry*ry + rz*rz);
  ggeo[gbase + Np*G01ID] = JW*(rx*sx + ry*sy + rz*sz);
  ggeo[gbase + Np*G02ID] = JW*(rx*tx + ry*ty + rz*tz);

  ggeo[gbase + Np*G11ID] = JW*(sx*sx + sy*sy + sz*sz);
  ggeo[gbase + Np*G12ID] = JW*(sx*tx + sy*ty + sz*tz);

  ggeo[gbase + Np*G22ID] = JW*(tx*tx + ty*ty + tz*tz);
  ggeo[gbase + Np*GWJID] = JW;

  // now do for cubvgeo
  // 1. interpolate Jacobian matrix to cubature nodes
  for(int m=0;m<cubNq;++m){
    for(int n=0;n<cubNq;++n){
      dfloat cIni = cubInterp[n*Nq+i];
      dfloat cImj = cubInterp[m*Nq+j];
      cxr[n+m*cubNq] += cIni*cImj*xr;
      cxs[n+m*cubNq] += cIni*cImj*xs;
      cyr[n+m*cubNq] += cIni*cImj*yr;
      cys[n+m*cubNq] += cIni*cImj*ys;
      czr[n+m*cubNq] += cIni*cImj*zr;
      czs[n+m*cubNq] += cIni*cImj*zs;
      cx[n+m*cubNq] += cIni*cImj*xij;
      cy[n+m*cubNq] += cIni*cImj*yij;
      cz[n+m*cubNq] += cIni*cImj*zij;
    }
  }
      }
    }


    for(int n=0;n<cubNq*cubNq;++n){

      dfloat rx = cys[n]*cz[n] - czs[n]*cy[n]; // dXds x X
      dfloat ry = czs[n]*cx[n] - cxs[n]*cz[n];
      dfloat rz = cxs[n]*cy[n] - cys[n]*cx[n];

      dfloat sx = czr[n]*cy[n] - cyr[n]*cz[n]; // -dXdr x X
      dfloat sy = cxr[n]*cz[n] - czr[n]*cx[n];
      dfloat sz = cyr[n]*cx[n] - cxr[n]*cy[n];

      dfloat tx = cyr[n]*czs[n] - czr[n]*cys[n]; // dXdr x dXds ~ X*|dXdr x dXds|/|X|
      dfloat ty = czr[n]*cxs[n] - cxr[n]*czs[n];
      dfloat tz = cxr[n]*cys[n] - cyr[n]*cxs[n];

      dfloat Gx = tx, Gy = ty, Gz = tz;

      dfloat J = cx[n]*tx + cy[n]*ty + cz[n]*tz;

      if(J<1e-8) {
        stringstream ss;
        ss << "Negative J found at element " << e << "\n";
        LIBP_ABORT(ss.str())
      }

      rx /= J;      sx /= J;      tx /= J;
      ry /= J;      sy /= J;      ty /= J;
      rz /= J;      sz /= J;      tz /= J;

      // use this for "volume" Jacobian
      J = sqrt(Gx*Gx+Gy*Gy+Gz*Gz);

      if(J<1e-8) {
        stringstream ss;
        ss << "Negative J found at element " << e << "\n";
        LIBP_ABORT(ss.str())
      }

      dfloat JW = J*cubw[n%cubNq]*cubw[n/cubNq];

      /* store geometric factors */
      int base = Nvgeo*cubNp*e + n;

      cubvgeo[base + cubNp*RXID] = rx;
      cubvgeo[base + cubNp*RYID] = ry;
      cubvgeo[base + cubNp*RZID] = rz;
      cubvgeo[base + cubNp*SXID] = sx;
      cubvgeo[base + cubNp*SYID] = sy;
      cubvgeo[base + cubNp*SZID] = sz;
      cubvgeo[base + cubNp*TXID] = tx;
      cubvgeo[base + cubNp*TYID] = ty;
      cubvgeo[base + cubNp*TZID] = tz;
      cubvgeo[base + cubNp*JID]  = J;
      cubvgeo[base + cubNp*JWID] = JW;
      cubvgeo[base + cubNp*IJWID] = 1./JW;
    }
  }

  HaloExchange(vgeo, Nvgeo*Np, ogsDfloat);
}
