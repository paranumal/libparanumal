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

void meshTri2D::PhysicalNodes(){

  x = (dfloat*) calloc((Nelements+totalHaloPairs)*Np,sizeof(dfloat));
  y = (dfloat*) calloc((Nelements+totalHaloPairs)*Np,sizeof(dfloat));
  z = (dfloat*) calloc((Nelements+totalHaloPairs)*Np,sizeof(dfloat)); // dummy

  dlong cnt = 0;
  for(dlong e=0;e<Nelements;++e){ /* for each element */

    dlong id = e*Nverts+0;

    dfloat xe1 = EX[id+0]; /* x-coordinates of vertices */
    dfloat xe2 = EX[id+1];
    dfloat xe3 = EX[id+2];

    dfloat ye1 = EY[id+0]; /* y-coordinates of vertices */
    dfloat ye2 = EY[id+1];
    dfloat ye3 = EY[id+2];

    for(int n=0;n<Np;++n){ /* for each node */

      /* (r,s) coordinates of interpolation nodes*/
      dfloat rn = r[n];
      dfloat sn = s[n];

      /* physical coordinate of interpolation node */
      x[cnt] = -0.5*(rn+sn)*xe1 + 0.5*(1+rn)*xe2 + 0.5*(1+sn)*xe3;
      y[cnt] = -0.5*(rn+sn)*ye1 + 0.5*(1+rn)*ye2 + 0.5*(1+sn)*ye3;
      ++cnt;
    }
  }

  HaloExchange(x, Np, ogsDfloat);
  HaloExchange(y, Np, ogsDfloat);

  // grab EX,EY,EZ from halo
  EX = (dfloat*) realloc(EX, (Nelements+totalHaloPairs)*Nverts*sizeof(dfloat));
  EY = (dfloat*) realloc(EY, (Nelements+totalHaloPairs)*Nverts*sizeof(dfloat));

  // send halo data and recv into extended part of arrays
  HaloExchange(EX, Nverts, ogsDfloat);
  HaloExchange(EY, Nverts, ogsDfloat);
}
