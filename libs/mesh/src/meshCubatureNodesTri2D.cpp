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

void meshTri2D::CubatureNodes(){

  cubx = (dfloat*) calloc(Nelements*cubNp,sizeof(dfloat));
  cuby = (dfloat*) calloc(Nelements*cubNp,sizeof(dfloat));

  dlong cnt = 0;
  for(dlong e=0;e<Nelements;++e){ /* for each element */

    dlong id = e*Nverts+0;

    dfloat xe1 = EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = EX[id+1];
    dfloat xe3 = EX[id+2];

    dfloat ye1 = EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = EY[id+1];
    dfloat ye3 = EY[id+2];

    for(int n=0;n<cubNp;++n){ /* for each node */

      /* (r,s) coordinates of interpolation nodes*/
      dfloat rn = cubr[n];
      dfloat sn = cubs[n];

      /* physical coordinate of interpolation node */
      cubx[cnt] = -0.5*(rn+sn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3;
      cuby[cnt] = -0.5*(rn+sn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3;
      ++cnt;
    }
  }

  o_cubx = device.malloc(Nelements*cubNp*sizeof(dfloat), cubx);
  o_cuby = device.malloc(Nelements*cubNp*sizeof(dfloat), cuby);
  o_cubz = o_cuby; // dummy to align with 3d

  //Face cubature
  intx = (dfloat*) calloc(Nelements*Nfaces*intNfp, sizeof(dfloat));
  inty = (dfloat*) calloc(Nelements*Nfaces*intNfp, sizeof(dfloat));
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      for(int n=0;n<intNfp;++n){
        dfloat ix = 0, iy = 0;
        for(int m=0;m<Nfp;++m){
          dlong vid = vmapM[m+f*Nfp+e*Nfp*Nfaces];
          dfloat xm = x[vid];
          dfloat ym = y[vid];
          dfloat Inm = intInterp[m+n*Nfp+f*intNfp*Nfp];
          ix += Inm*xm;
          iy += Inm*ym;
        }
        dlong id = n + f*intNfp + e*Nfaces*intNfp;
        intx[id] = ix;
        inty[id] = iy;
      }
    }
  }

  o_intx = device.malloc(Nelements*Nfaces*intNfp*sizeof(dfloat), intx);
  o_inty = device.malloc(Nelements*Nfaces*intNfp*sizeof(dfloat), inty);
  o_intz = o_inty; // dummy to align with 3d
}
