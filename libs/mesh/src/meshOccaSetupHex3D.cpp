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

void meshHex3D::OccaSetup(){

  this->mesh3D::OccaSetup();

  //lumped mass matrix
  MM = (dfloat *) calloc(Np*Np, sizeof(dfloat));
  for (int k=0;k<Nq;k++) {
    for (int j=0;j<Nq;j++) {
      for (int i=0;i<Nq;i++) {
        int n = i+j*Nq+k*Nq*Nq;
        MM[n+n*Np] = gllw[i]*gllw[j]*gllw[k];
      }
    }
  }

  //build inverse of mass matrix
  invMM = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  for (int n=0;n<Np*Np;n++)
    invMM[n] = MM[n];
  matrixInverse(Np,invMM);

  // build trilinear geometric factors for hexes
  if(settings.compareSetting("ELEMENT MAP", "AFFINE")){
    // pack gllz, gllw, and elementwise EXYZ
    hlong Nxyz = Nelements*dim*Nverts;
    EXYZ  = (dfloat*) calloc(Nxyz, sizeof(dfloat));
    gllzw = (dfloat*) calloc(2*Nq, sizeof(dfloat));

    int sk = 0;
    for(int n=0;n<Nq;++n)
      gllzw[sk++] = gllz[n];
    for(int n=0;n<Nq;++n)
      gllzw[sk++] = gllw[n];

    sk = 0;
    for(hlong e=0;e<Nelements;++e){
      for(int v=0;v<Nverts;++v)
        EXYZ[sk++] = EX[e*Nverts+v];
      for(int v=0;v<Nverts;++v)
        EXYZ[sk++] = EY[e*Nverts+v];
      for(int v=0;v<Nverts;++v)
        EXYZ[sk++] = EZ[e*Nverts+v];
    }

    // nodewise ggeo with element coordinates and gauss node info
    o_EXYZ  = device.malloc(Nxyz*sizeof(dfloat), EXYZ);
    o_gllzw = device.malloc(2*Nq*sizeof(dfloat), gllzw);
  }

  ggeoNoJW = (dfloat*) calloc(Np*Nelements*6,sizeof(dfloat));
  for(int e=0;e<Nelements;++e){
    for(int n=0;n<Np;++n){
      ggeoNoJW[e*Np*6 + n + 0*Np] = ggeo[e*Np*Nggeo + n + G00ID*Np];
      ggeoNoJW[e*Np*6 + n + 1*Np] = ggeo[e*Np*Nggeo + n + G01ID*Np];
      ggeoNoJW[e*Np*6 + n + 2*Np] = ggeo[e*Np*Nggeo + n + G02ID*Np];
      ggeoNoJW[e*Np*6 + n + 3*Np] = ggeo[e*Np*Nggeo + n + G11ID*Np];
      ggeoNoJW[e*Np*6 + n + 4*Np] = ggeo[e*Np*Nggeo + n + G12ID*Np];
      ggeoNoJW[e*Np*6 + n + 5*Np] = ggeo[e*Np*Nggeo + n + G22ID*Np];
    }
  }
  o_ggeoNoJW = device.malloc(Np*Nelements*6*sizeof(dfloat), ggeoNoJW);


  o_MM = device.malloc(Np*Np*sizeof(dfloat), MM);

  o_D = device.malloc(Nq*Nq*sizeof(dfloat), D);

  o_Dmatrices = device.malloc(Nq*Nq*sizeof(dfloat), D);
  o_Smatrices = device.malloc(Nq*Nq*sizeof(dfloat), D); //dummy

  o_sMT = device.malloc(1*sizeof(dfloat)); //dummy
  o_LIFTT = device.malloc(1*sizeof(dfloat)); // dummy

  o_vgeo =
    device.malloc((Nelements+totalHaloPairs)*Np*Nvgeo*sizeof(dfloat),
                        vgeo);

  o_sgeo =
    device.malloc(Nelements*Nfaces*Nfp*Nsgeo*sizeof(dfloat),
                        sgeo);

  o_ggeo =
    device.malloc(Nelements*Np*Nggeo*sizeof(dfloat),
      ggeo);
}
