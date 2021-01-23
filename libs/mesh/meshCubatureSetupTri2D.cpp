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
#include "mesh/mesh3D.hpp"

void meshTri3D::CubatureSetup(){

  CubatureSetup(N, " ");
}

void meshTri3D::CubatureSetup(int _cubN, const char *_cubatureType){
  mesh_t *mesh_p = (mesh_t*) this;
  meshTri2D* trimesh = (meshTri2D*) mesh_p;
  trimesh->meshTri2D::CubatureSetup(_cubN, _cubatureType);
}

void meshTri2D::CubatureSetup(){

  CubatureSetup(N, " ");
}

void meshTri2D::CubatureSetup(int _cubN, const char *_cubatureType){

  cubatureType = strdup(_cubatureType);

  /* Cubature data */
  cubN = 2*_cubN; //cubature order
  CubatureNodesTri2D(cubN, &cubNp, &cubr, &cubs, &cubw);

  cubInterp = (dfloat *) malloc(Np*cubNp*sizeof(dfloat));
  InterpolationMatrixTri2D(N, Np, r, s, cubNp, cubr, cubs, cubInterp);

  //cubature project cubProject = M^{-1} * cubInterp^T
  // Defined such that cubProject * cubW * cubInterp = Identity
  cubProject = (dfloat*) calloc(cubNp*Np, sizeof(dfloat));
  CubaturePmatrixTri2D(N, Np, r, s, cubNp, cubr, cubs, cubProject);

  //cubature derivates matrices, cubD: differentiate on cubature nodes
  // we dont use cubD on Tris/Tets  so skip computing

  // Instead, it's cheaper to:
  // make weak cubature derivatives cubPDT = cubProject * cubD^T
  cubPDT  = (dfloat*) calloc(2*cubNp*Np, sizeof(dfloat));
  cubPDrT = cubPDT + 0*cubNp*Np;
  cubPDsT = cubPDT + 1*cubNp*Np;
  CubatureWeakDmatricesTri2D(N, Np, r, s, cubNp, cubr, cubs, cubPDrT, cubPDsT);

  // cubN+1 point Gauss-Legendre quadrature for surface integrals
  cubNq  = cubN+1;
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

  // build transposes (we hold matrices as column major on device)
  dfloat *cubProjectT = (dfloat*) calloc(cubNp*Np, sizeof(dfloat));
  dfloat *cubInterpT   = (dfloat*) calloc(cubNp*Np, sizeof(dfloat));
  matrixTranspose(cubNp, Np, cubInterp, Np, cubInterpT, cubNp);
  matrixTranspose(Np, cubNp, cubProject, cubNp, cubProjectT, Np);

  //pre-multiply cubProject by W on device
  for(int n=0;n<cubNp;++n){
    for(int m=0;m<Np;++m){
      cubProjectT[m+n*Np] *= cubw[n];
    }
  }

  dfloat *cubPDTT = (dfloat*) calloc(2*cubNp*Np, sizeof(dfloat));
  dfloat *cubPDrTT = cubPDTT + 0*cubNp*Np;
  dfloat *cubPDsTT = cubPDTT + 1*cubNp*Np;
  matrixTranspose(Np, cubNp, cubPDrT, cubNp, cubPDrTT, Np);
  matrixTranspose(Np, cubNp, cubPDsT, cubNp, cubPDsTT, Np);

  //pre-multiply cubPDT by W on device
  for(int n=0;n<cubNp;++n){
    for(int m=0;m<Np;++m){
      cubPDrTT[m+n*Np] *= cubw[n];
      cubPDsTT[m+n*Np] *= cubw[n];
    }
  }

  // build surface integration matrix transposes
  dfloat *intLIFTT = (dfloat*) calloc(Np*Nfaces*intNfp, sizeof(dfloat));
  dfloat *intInterpT = (dfloat*) calloc(Nfp*Nfaces*intNfp, sizeof(dfloat));
  matrixTranspose(Np, Nfaces*intNfp, intLIFT, Nfaces*intNfp, intLIFTT, Np);
  matrixTranspose(Nfaces*intNfp, Nfp, intInterp, Nfp, intInterpT, Nfaces*intNfp);

  o_cubvgeo = o_vgeo;// dummy
  o_cubsgeo = o_sgeo; //dummy cubature geo factors

  o_cubInterp   = platform.malloc(Np*cubNp*sizeof(dfloat), cubInterpT);
  o_cubProject = platform.malloc(Np*cubNp*sizeof(dfloat), cubProjectT);

  o_cubPDT = platform.malloc(2*cubNp*Np*sizeof(dfloat), cubPDTT);
  o_cubD = o_cubPDT; //dummy

  o_intInterp = platform.malloc(Nfp*Nfaces*intNfp*sizeof(dfloat), intInterpT);
  o_intLIFT = platform.malloc(Np*Nfaces*intNfp*sizeof(dfloat), intLIFTT);

  free(cubPDTT);
  free(cubProjectT);
  free(cubInterpT);
  free(intLIFTT);
  free(intInterpT);
}
