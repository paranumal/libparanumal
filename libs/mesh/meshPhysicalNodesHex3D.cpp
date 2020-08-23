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

void meshHex3D::PhysicalNodes(){

  x = (dfloat*) calloc((Nelements+totalHaloPairs)*Np,sizeof(dfloat));
  y = (dfloat*) calloc((Nelements+totalHaloPairs)*Np,sizeof(dfloat));
  z = (dfloat*) calloc((Nelements+totalHaloPairs)*Np,sizeof(dfloat));

  dlong cnt = 0;
  for(dlong e=0;e<Nelements;++e){ /* for each element */

    dlong id = e*Nverts;

    dfloat xe1 = EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = EX[id+1];
    dfloat xe3 = EX[id+2];
    dfloat xe4 = EX[id+3];
    dfloat xe5 = EX[id+4];
    dfloat xe6 = EX[id+5];
    dfloat xe7 = EX[id+6];
    dfloat xe8 = EX[id+7];

    dfloat ye1 = EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = EY[id+1];
    dfloat ye3 = EY[id+2];
    dfloat ye4 = EY[id+3];
    dfloat ye5 = EY[id+4];
    dfloat ye6 = EY[id+5];
    dfloat ye7 = EY[id+6];
    dfloat ye8 = EY[id+7];

    dfloat ze1 = EZ[id+0]; /* z-coordinates of vertices */
    dfloat ze2 = EZ[id+1];
    dfloat ze3 = EZ[id+2];
    dfloat ze4 = EZ[id+3];
    dfloat ze5 = EZ[id+4];
    dfloat ze6 = EZ[id+5];
    dfloat ze7 = EZ[id+6];
    dfloat ze8 = EZ[id+7];

    for(int n=0;n<Np;++n){ /* for each node */

      /* (r,s,t) coordinates of interpolation nodes*/
      dfloat rn = r[n];
      dfloat sn = s[n];
      dfloat tn = t[n];

      /* physical coordinate of interpolation node */
      x[cnt] =
        +0.125*(1-rn)*(1-sn)*(1-tn)*xe1
        +0.125*(1+rn)*(1-sn)*(1-tn)*xe2
        +0.125*(1+rn)*(1+sn)*(1-tn)*xe3
        +0.125*(1-rn)*(1+sn)*(1-tn)*xe4
        +0.125*(1-rn)*(1-sn)*(1+tn)*xe5
        +0.125*(1+rn)*(1-sn)*(1+tn)*xe6
        +0.125*(1+rn)*(1+sn)*(1+tn)*xe7
        +0.125*(1-rn)*(1+sn)*(1+tn)*xe8;

      y[cnt] =
        +0.125*(1-rn)*(1-sn)*(1-tn)*ye1
        +0.125*(1+rn)*(1-sn)*(1-tn)*ye2
        +0.125*(1+rn)*(1+sn)*(1-tn)*ye3
        +0.125*(1-rn)*(1+sn)*(1-tn)*ye4
        +0.125*(1-rn)*(1-sn)*(1+tn)*ye5
        +0.125*(1+rn)*(1-sn)*(1+tn)*ye6
        +0.125*(1+rn)*(1+sn)*(1+tn)*ye7
        +0.125*(1-rn)*(1+sn)*(1+tn)*ye8;

      z[cnt] =
        +0.125*(1-rn)*(1-sn)*(1-tn)*ze1
        +0.125*(1+rn)*(1-sn)*(1-tn)*ze2
        +0.125*(1+rn)*(1+sn)*(1-tn)*ze3
        +0.125*(1-rn)*(1+sn)*(1-tn)*ze4
        +0.125*(1-rn)*(1-sn)*(1+tn)*ze5
        +0.125*(1+rn)*(1-sn)*(1+tn)*ze6
        +0.125*(1+rn)*(1+sn)*(1+tn)*ze7
        +0.125*(1-rn)*(1+sn)*(1+tn)*ze8;

      ++cnt;
    }
  }

  halo->Exchange(x, Np, ogs_dfloat);
  halo->Exchange(y, Np, ogs_dfloat);
  halo->Exchange(z, Np, ogs_dfloat);
}
