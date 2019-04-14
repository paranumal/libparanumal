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

void meshQuad2D::GeometricFactors(){

  /* unified storage array for geometric factors */
  Nvgeo = 7;

  /* note that we have volume geometric factors for each node */
  vgeo = (dfloat*) calloc(Nelements*Nvgeo*Np, sizeof(dfloat));

  cubvgeo = (dfloat*) calloc(Nelements*Nvgeo*cubNp, sizeof(dfloat));

  /* number of second order geometric factors */
  Nggeo = 4;
  ggeo = (dfloat*) calloc(Nelements*Nggeo*Np, sizeof(dfloat));

  for(dlong e=0;e<Nelements;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    dlong id = e*Nverts;

    dfloat *xe = EX + id;
    dfloat *ye = EY + id;

    for(int n=0;n<Np;++n){

      /* local node coordinates */
      dfloat rn = r[n];
      dfloat sn = s[n];

      /* Jacobian matrix */
      dfloat xr = 0.25*( (1-sn)*(xe[1]-xe[0]) + (1+sn)*(xe[2]-xe[3]) );
      dfloat xs = 0.25*( (1-rn)*(xe[3]-xe[0]) + (1+rn)*(xe[2]-xe[1]) );
      dfloat yr = 0.25*( (1-sn)*(ye[1]-ye[0]) + (1+sn)*(ye[2]-ye[3]) );
      dfloat ys = 0.25*( (1-rn)*(ye[3]-ye[0]) + (1+rn)*(ye[2]-ye[1]) );

      /* compute geometric factors for affine coordinate transform*/
      dfloat J = xr*ys - xs*yr;

      if(J<1e-8) {
        stringstream ss;
        ss << "Negative J found at element " << e << "\n";
        LIBP_ABORT(ss.str())
      }
      dfloat rx =  ys/J;
      dfloat ry = -xs/J;
      dfloat sx = -yr/J;
      dfloat sy =  xr/J;

      int i = n%Nq;
      int j = n/Nq;
      dfloat JW = J*gllw[i]*gllw[j];

      /* store geometric factors */
      vgeo[Nvgeo*Np*e + n + Np*RXID] = rx;
      vgeo[Nvgeo*Np*e + n + Np*RYID] = ry;
      vgeo[Nvgeo*Np*e + n + Np*SXID] = sx;
      vgeo[Nvgeo*Np*e + n + Np*SYID] = sy;
      vgeo[Nvgeo*Np*e + n + Np*JID]  = J;
      vgeo[Nvgeo*Np*e + n + Np*JWID] = JW;
      vgeo[Nvgeo*Np*e + n + Np*IJWID] = 1./JW;

      /* store second order geometric factors */
      ggeo[Nggeo*Np*e + n + Np*G00ID] = JW*(rx*rx + ry*ry);
      ggeo[Nggeo*Np*e + n + Np*G01ID] = JW*(rx*sx + ry*sy);
      ggeo[Nggeo*Np*e + n + Np*G11ID] = JW*(sx*sx + sy*sy);
      ggeo[Nggeo*Np*e + n + Np*GWJID] = JW;
    }

    //geometric data for quadrature
    for(int j=0;j<cubNq;++j){
      for(int i=0;i<cubNq;++i){

        dfloat rn = cubr[i];
        dfloat sn = cubr[j];

        /* Jacobian matrix */
        dfloat xr = 0.25*( (1-sn)*(xe[1]-xe[0]) + (1+sn)*(xe[2]-xe[3]) );
        dfloat xs = 0.25*( (1-rn)*(xe[3]-xe[0]) + (1+rn)*(xe[2]-xe[1]) );
        dfloat yr = 0.25*( (1-sn)*(ye[1]-ye[0]) + (1+sn)*(ye[2]-ye[3]) );
        dfloat ys = 0.25*( (1-rn)*(ye[3]-ye[0]) + (1+rn)*(ye[2]-ye[1]) );

        /* compute geometric factors for affine coordinate transform*/
        dfloat J = xr*ys - xs*yr;

        if(J<1e-8) {
          stringstream ss;
          ss << "Negative J found at element " << e << "\n";
          LIBP_ABORT(ss.str())
        }
        dfloat rx =  ys/J;
        dfloat ry = -xs/J;
        dfloat sx = -yr/J;
        dfloat sy =  xr/J;

        dfloat JW = J*cubw[i]*cubw[j];

        /* store geometric factors */
        dlong base = Nvgeo*cubNp*e + i + j*cubNq;
        cubvgeo[base + cubNp*RXID] = rx;
        cubvgeo[base + cubNp*RYID] = ry;
        cubvgeo[base + cubNp*SXID] = sx;
        cubvgeo[base + cubNp*SYID] = sy;
        cubvgeo[base + cubNp*JID]  = J;
        cubvgeo[base + cubNp*JWID] = JW;
        cubvgeo[base + cubNp*IJWID] = 1./JW;
      }
    }
  }
}
