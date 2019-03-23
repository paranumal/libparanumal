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

// custom geometric factors specialized for 3D tri on sphere

void meshTri3D::GeometricFactors(){

  /* unified storage array for geometric factors */
  Nvgeo = 12; //

  /* note that we have volume geometric factors for each node */
  vgeo = (dfloat*) calloc(Nelements*Nvgeo*Np, sizeof(dfloat));

  /* number of second order geometric factors */
  Nggeo = 7;
  ggeo = (dfloat*) calloc(Nelements*Nggeo, sizeof(dfloat));

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

      if(J<1e-8) {
        stringstream ss;
        ss << "Negative J found at element " << e << "\n";
        LIBP_ABORT(ss.str())
      }

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

      if(J<1e-8) {
        stringstream ss;
        ss << "Negative J found at element " << e << "\n";
        LIBP_ABORT(ss.str())
      }

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

      ggeo[Nggeo*e + G00ID] = J*(rx*rx + ry*ry + rz*rz);
      ggeo[Nggeo*e + G01ID] = J*(rx*sx + ry*sy + rz*sz);
      ggeo[Nggeo*e + G11ID] = J*(sx*sx + sy*sy + sz*sz);
      ggeo[Nggeo*e + GWJID]  = J;

    }
  }
}
