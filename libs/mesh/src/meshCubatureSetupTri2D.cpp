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
#include "mesh3D.hpp"

void meshTri3D::CubatureSetup(){
  mesh_t *mesh_p = (mesh_t*) this;
  meshTri2D* trimesh = (meshTri2D*) mesh_p;
  trimesh->meshTri2D::CubatureSetup();
}

void meshTri2D::CubatureSetup(){

  /* Cubature data */
  cubN = 2*N; //cubature order
  CubatureNodesTri2D(cubN, &cubNp, &cubr, &cubs, &cubw);

  cubInterp = (dfloat *) malloc(Np*cubNp*sizeof(dfloat));
  InterpolationMatrixTri2D(N, Np, r, s, cubNp, cubr, cubs, cubInterp);

  cubDrW = (dfloat*) calloc(cubNp*Np, sizeof(dfloat));
  cubDsW = (dfloat*) calloc(cubNp*Np, sizeof(dfloat));
  cubProject = (dfloat*) calloc(cubNp*Np, sizeof(dfloat));
  CubatureWeakDmatricesTri2D(N, Np, r, s, cubNp, cubr, cubs, cubw,
                             cubDrW, cubDsW, cubProject);

  // cubN+1 point Gauss-Legendre quadrature for surface integrals
  cubNq = cubN+1;
  cubNfp = cubN+1;
  intNfp = cubN+1;
  intr = (dfloat *) malloc(cubNfp*sizeof(dfloat));
  intw = (dfloat *) malloc(cubNfp*sizeof(dfloat));
  JacobiGQ(0, 0, cubN, intr, intw);

  intInterp = (dfloat*) calloc(intNfp*Nfaces*Nfp, sizeof(dfloat));
  intLIFT = (dfloat*) calloc(Nfaces*intNfp*Np, sizeof(dfloat));
  CubatureSurfaceMatricesTri2D(N, Np, r, s, faceNodes, intNfp, intr, intw,
                               intInterp, intLIFT);

  // add compile time constants to kernels
  props["defines/" "p_cubNq"]= cubNq;
  props["defines/" "p_cubNp"]= cubNp;
  props["defines/" "p_intNfp"]= intNfp;
  props["defines/" "p_intNfpNfaces"]= intNfp*Nfaces;
  props["defines/" "p_cubNfp"]= cubNfp;

  // build volume cubature matrix transposes
  int cubNpBlocked = Np*((cubNp+Np-1)/Np);
  dfloat *cubDrWT = (dfloat*) calloc(cubNpBlocked*Np, sizeof(dfloat));
  dfloat *cubDsWT = (dfloat*) calloc(cubNpBlocked*Np, sizeof(dfloat));
  dfloat *cubDrsWT = (dfloat*) calloc(2*cubNp*Np, sizeof(dfloat));
  dfloat *cubProjectT = (dfloat*) calloc(cubNp*Np, sizeof(dfloat));
  dfloat *cubInterpT = (dfloat*) calloc(cubNp*Np, sizeof(dfloat));
  for(int n=0;n<Np;++n){
    for(int m=0;m<cubNp;++m){
      cubDrWT[n+m*Np] = cubDrW[n*cubNp+m];
      cubDsWT[n+m*Np] = cubDsW[n*cubNp+m];

      cubDrsWT[n+m*Np] = cubDrW[n*cubNp+m];
      cubDrsWT[n+m*Np+cubNp*Np] = cubDsW[n*cubNp+m];

      cubProjectT[n+m*Np] = cubProject[n*cubNp+m];
      cubInterpT[m+n*cubNp] = cubInterp[m*Np+n];
    }
  }

  // build surface integration matrix transposes
  dfloat *intLIFTT = (dfloat*) calloc(Np*Nfaces*intNfp, sizeof(dfloat));
  dfloat *intInterpT = (dfloat*) calloc(Nfp*Nfaces*intNfp, sizeof(dfloat));
  for(int n=0;n<Np;++n){
    for(int m=0;m<Nfaces*intNfp;++m){
      intLIFTT[n+m*Np] = intLIFT[n*intNfp*Nfaces+m];
    }
  }
  for(int n=0;n<intNfp*Nfaces;++n){
    for(int m=0;m<Nfp;++m){
      intInterpT[n+m*Nfaces*intNfp] = intInterp[n*Nfp + m];
    }
  }

  o_cubvgeo = o_vgeo;// dummy

  o_cubsgeo = o_sgeo; //dummy cubature geo factors

  o_cubInterpT =
    device.malloc(Np*cubNp*sizeof(dfloat),
        cubInterpT);

  o_cubProjectT =
    device.malloc(Np*cubNp*sizeof(dfloat),
        cubProjectT);

  o_cubDrWT =
    device.malloc(Np*cubNp*sizeof(dfloat),
        cubDrWT);

  o_cubDsWT =
    device.malloc(Np*cubNp*sizeof(dfloat),
        cubDsWT);

  o_cubDtWT =
    device.malloc(Np*cubNp*sizeof(dfloat),
		  cubDsWT); // dummy to align with 3d

  o_cubDWmatrices = device.malloc(2*cubNp*Np*sizeof(dfloat), cubDrsWT);

  o_intInterpT =
    device.malloc(Nfp*Nfaces*intNfp*sizeof(dfloat),
        intInterpT);

  o_intLIFTT =
    device.malloc(Np*Nfaces*intNfp*sizeof(dfloat),
        intLIFTT);

  free(cubDrWT);
  free(cubDsWT);
  free(cubDrsWT);
  free(cubProjectT);
  free(cubInterpT);
}
