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

void meshTri2D::GeometricFactors(){

  /* unified storage array for geometric factors */
  Nvgeo = 5;
  vgeo = (dfloat*) calloc((Nelements+totalHaloPairs)*Nvgeo, sizeof(dfloat));

  /* number of second order geometric factors */
  Nggeo = 4;
  ggeo = (dfloat*) calloc(Nelements*Nggeo, sizeof(dfloat));


  for(dlong e=0;e<Nelements;++e){ /* for each element */

    /* find vertex indices and physical coordinates */
    dlong id = e*Nverts+0;

    dfloat xe1 = EX[id+0];
    dfloat xe2 = EX[id+1];
    dfloat xe3 = EX[id+2];

    dfloat ye1 = EY[id+0];
    dfloat ye2 = EY[id+1];
    dfloat ye3 = EY[id+2];

    /* compute geometric factors for affine coordinate transform*/
    dfloat J = 0.25*((xe2-xe1)*(ye3-ye1) - (xe3-xe1)*(ye2-ye1));

    if(J<0) {
      stringstream ss;
      ss << "Negative J found at element " << e << "\n";
      LIBP_ABORT(ss.str())
    }
    dfloat rx =  (0.5/J)*(ye3-ye1);
    dfloat ry = -(0.5/J)*(xe3-xe1);
    dfloat sx = -(0.5/J)*(ye2-ye1);
    dfloat sy =  (0.5/J)*(xe2-xe1);

    /* store geometric factors */
    vgeo[Nvgeo*e + RXID] = rx;
    vgeo[Nvgeo*e + RYID] = ry;
    vgeo[Nvgeo*e + SXID] = sx;
    vgeo[Nvgeo*e + SYID] = sy;
    vgeo[Nvgeo*e +  JID] = J;

    /* store second order geometric factors */
    ggeo[Nggeo*e + G00ID] = J*(rx*rx + ry*ry);
    ggeo[Nggeo*e + G01ID] = J*(rx*sx + ry*sy);
    ggeo[Nggeo*e + G11ID] = J*(sx*sx + sy*sy);
    ggeo[Nggeo*e + GWJID]  = J;
  }

  halo->Exchange(vgeo, Nvgeo, ogs_dfloat);
}
