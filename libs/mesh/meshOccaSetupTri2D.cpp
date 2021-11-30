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
#include "mesh/mesh2D.hpp"

void meshTri2D::OccaSetup(){

  this->mesh2D::OccaSetup();

  // build transposes (we hold matrices as column major on device)
  dfloat *DT = (dfloat*) calloc(2*Np*Np, sizeof(dfloat));
  dfloat *DrT = DT + 0*Np*Np;
  dfloat *DsT = DT + 1*Np*Np;
  matrixTranspose(Np, Np, Dr, Np, DrT, Np);
  matrixTranspose(Np, Np, Ds, Np, DsT, Np);

  // build transposes (we hold matrices as column major on device)
  dfloat *DWT = (dfloat*) calloc(2*Np*Np, sizeof(dfloat));
  dfloat *DWrT = DWT + 0*Np*Np;
  dfloat *DWsT = DWT + 1*Np*Np;
  matrixTranspose(Np, Np, DWr, Np, DWrT, Np);
  matrixTranspose(Np, Np, DWs, Np, DWsT, Np);

  dfloat *LIFTT = (dfloat*) calloc(Np*Nfaces*Nfp, sizeof(dfloat));
  matrixTranspose(Np, Nfp*Nfaces, LIFT, Nfp*Nfaces, LIFTT, Np);

  dfloat *sMT = (dfloat *) calloc(Np*Nfaces*Nfp,sizeof(dfloat));
  matrixTranspose(Np, Nfp*Nfaces, sM, Nfp*Nfaces, sMT, Np);

  dfloat *ST = (dfloat*) calloc(3*Np*Np, sizeof(dfloat));
  dfloat *SrrT = ST + 0*Np*Np;
  dfloat *SrsT = ST + 1*Np*Np;
  dfloat *SssT = ST + 2*Np*Np;
  matrixTranspose(Np, Np, Srr, Np, SrrT, Np);
  matrixTranspose(Np, Np, Srs, Np, SrsT, Np);
  matrixTranspose(Np, Np, Sss, Np, SssT, Np);

  o_D  = platform.malloc(2*Np*Np*sizeof(dfloat), DT);
  o_DW = platform.malloc(2*Np*Np*sizeof(dfloat), DWT);
  o_MM = platform.malloc(Np*Np*sizeof(dfloat), MM); //MM is symmetric

  o_sM = platform.malloc(Np*Nfaces*Nfp*sizeof(dfloat), sMT);

  o_LIFT = platform.malloc(Np*Nfaces*Nfp*sizeof(dfloat), LIFTT);

  o_S = platform.malloc(3*Np*Np*sizeof(dfloat), ST);

  o_vgeo = platform.malloc((Nelements+totalHaloPairs)*Nvgeo*sizeof(dfloat), vgeo);
  o_sgeo = platform.malloc(Nelements*Nfaces*Nsgeo*sizeof(dfloat), sgeo);
  o_ggeo = platform.malloc(Nelements*Nggeo*sizeof(dfloat), ggeo);

  free(DT);
  free(DWT);
  free(LIFTT);
  free(sMT);
  free(ST);
}
