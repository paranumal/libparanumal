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

void meshQuad2D::OccaSetup(){

  this->mesh2D::OccaSetup();

  //lumped mass matrix
  MM = (dfloat *) calloc(Np*Np, sizeof(dfloat));
  for (int j=0;j<Nq;j++) {
    for (int i=0;i<Nq;i++) {
      int n = i+j*Nq;
      MM[n+n*Np] = gllw[i]*gllw[j];
    }
  }

  //build inverse of mass matrix
  invMM = (dfloat *) calloc(Np*Np,sizeof(dfloat));
  for (int n=0;n<Np*Np;n++)
    invMM[n] = MM[n];
  matrixInverse(Np,invMM);

  dfloat *cubDWT = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  dfloat *cubProjectT = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  dfloat *cubInterpT = (dfloat*) calloc(cubNq*Nq, sizeof(dfloat));
  for(int n=0;n<Nq;++n){
    for(int m=0;m<cubNq;++m){
      cubDWT[n+m*Nq] = cubDW[n*cubNq+m];
      cubProjectT[n+m*Nq] = cubProject[n*cubNq+m];
      cubInterpT[m+n*cubNq] = cubInterp[m*Nq+n];
    }
  }

  intx = (dfloat*) calloc(Nelements*Nfaces*cubNq, sizeof(dfloat));
  inty = (dfloat*) calloc(Nelements*Nfaces*cubNq, sizeof(dfloat));
  for(dlong e=0;e<Nelements;++e){
    for(int f=0;f<Nfaces;++f){
      for(int n=0;n<cubNq;++n){
        dfloat ix = 0, iy = 0;
        for(int m=0;m<Nq;++m){
          dlong vid = vmapM[m+f*Nfp+e*Nfp*Nfaces];
          dfloat xm = x[vid];
          dfloat ym = y[vid];

          dfloat Inm = cubInterp[m+n*Nq];
          ix += Inm*xm;
          iy += Inm*ym;
        }
        dlong id = n + f*cubNq + e*Nfaces*cubNq;
        intx[id] = ix;
        inty[id] = iy;
      }
    }
  }

  o_D = device.malloc(Nq*Nq*sizeof(dfloat), D);

  o_Dmatrices = device.malloc(Nq*Nq*sizeof(dfloat), D);
  o_Smatrices = device.malloc(Nq*Nq*sizeof(dfloat), D); //dummy

  o_MM = device.malloc(Np*Np*sizeof(dfloat), MM);

  o_sMT = device.malloc(1*sizeof(dfloat)); //dummy

  o_vgeo =
    device.malloc((Nelements+totalHaloPairs)*Nvgeo*Np*sizeof(dfloat),
        vgeo);
  o_sgeo =
    device.malloc(Nelements*Nfaces*Nfp*Nsgeo*sizeof(dfloat),
        sgeo);
  o_ggeo =
    device.malloc(Nelements*Np*Nggeo*sizeof(dfloat),
        ggeo);

  o_cubvgeo =
    device.malloc(Nelements*Nvgeo*cubNp*sizeof(dfloat),
        cubvgeo);

  o_cubsgeo =
    device.malloc(Nelements*Nfaces*cubNq*Nsgeo*sizeof(dfloat),
        cubsgeo);

  o_cubInterpT =
    device.malloc(Nq*cubNq*sizeof(dfloat),
        cubInterpT);

  o_cubProjectT =
    device.malloc(Nq*cubNq*sizeof(dfloat),
        cubProjectT);

  o_cubDWT =
    device.malloc(Nq*cubNq*sizeof(dfloat),
        cubDWT);

  o_cubDWmatrices = device.malloc(cubNq*Nq*sizeof(dfloat), cubDWT);

  o_LIFTT =
    device.malloc(1*sizeof(dfloat)); // dummy

  o_intx =
    device.malloc(Nelements*Nfaces*cubNq*sizeof(dfloat),
        intx);

  o_inty =
    device.malloc(Nelements*Nfaces*cubNq*sizeof(dfloat),
        inty);

  o_intz =
    device.malloc(Nelements*Nfaces*cubNq*sizeof(dfloat),
        inty);

  //dummy quadrature lifter operators
  o_intInterpT = device.malloc(cubNq*Nq*sizeof(dfloat));
  o_intInterpT.copyFrom(o_cubInterpT);

  o_intLIFTT = device.malloc(cubNq*Nq*sizeof(dfloat));
  o_intLIFTT.copyFrom(o_cubProjectT);

}
