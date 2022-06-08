/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

namespace libp {

void mesh_t::SetupPmlBoxTri2D(){

  // find a factorization size = size_x*size_y such that
  //  size_x>=size_y and are 'close' to one another
  int size_x, size_y;
  Factor2(size, size_x, size_y);

  //determine (x,y) rank coordinates for this processes
  int rank_x=-1, rank_y=-1;
  RankDecomp2(size_x, size_y,
              rank_x, rank_y,
              rank);

  //get global size from settings
  dlong NX, NY;
  settings.getSetting("BOX GLOBAL NX", NX);
  settings.getSetting("BOX GLOBAL NY", NY);

  //number of local elements in each dimension
  dlong nx, ny;
  settings.getSetting("BOX NX", nx);
  settings.getSetting("BOX NY", ny);

  if (NX*NY <= 0) { //if the user hasn't given global sizes
    //set global size by multiplying local size by grid dims
    NX = nx * size_x;
    NY = ny * size_y;
    settings.changeSetting("BOX GLOBAL NX", std::to_string(NX));
    settings.changeSetting("BOX GLOBAL NY", std::to_string(NY));
  } else {
    //WARNING setting global sizes on input overrides any local sizes provided
    nx = NX/size_x + ((rank_x < (NX % size_x)) ? 1 : 0);
    ny = NY/size_y + ((rank_y < (NY % size_y)) ? 1 : 0);
  }

  int boundaryFlag;
  settings.getSetting("BOX BOUNDARY FLAG", boundaryFlag);

  const int periodicFlag = (boundaryFlag == -1) ? 1 : 0;
  LIBP_ABORT("Periodic boundary unsupported for PMLBOX mesh.",
             periodicFlag);

  //local grid physical sizes
  dfloat DIMX, DIMY;
  settings.getSetting("BOX DIMX", DIMX);
  settings.getSetting("BOX DIMY", DIMY);

  //relative PML width
  const dfloat pmlScale = 0.2;

  //element sizes
  dfloat dx = DIMX/NX;
  dfloat dy = DIMY/NY;

  dlong offset_x = rank_x*(NX/size_x) + std::min(rank_x, (NX % size_x));
  dlong offset_y = rank_y*(NY/size_y) + std::min(rank_y, (NY % size_y));

  //local grid physical sizes
  dfloat dimx = nx*dx;
  dfloat dimy = ny*dy;

  //pml width
  dfloat pmlWidthx = pmlScale*DIMX;
  dfloat pmlWidthy = pmlScale*DIMY;

  //bottom corner of physical domain
  dfloat X0 = -DIMX/2.0 + offset_x*dx;
  dfloat Y0 = -DIMY/2.0 + offset_y*dy;

  //try and have roughly the same resolution in the pml region
  dlong pmlNx, pmlNy;
  pmlNx = ceil(pmlWidthx/dx);
  pmlNy = ceil(pmlWidthy/dy);

  dfloat pmldx, pmldy;
  pmldx = pmlWidthx/pmlNx;
  pmldy = pmlWidthy/pmlNy;

  NX += 2*pmlNx;
  NY += 2*pmlNy;

  //global number of nodes in each dimension
  hlong NnX = periodicFlag ? NX : NX+1; //lose a node when periodic (repeated node)
  hlong NnY = periodicFlag ? NY : NY+1; //lose a node when periodic (repeated node)

  // build an nx x ny x nz box grid
  Nnodes = NnX*NnY; //global node count

  Nelements = 2*nx*ny; //local

  //pml sides
  if (rank_x==0)        Nelements+=2*ny*pmlNx;
  if (rank_x==size_x-1) Nelements+=2*ny*pmlNx;
  if (rank_y==0)        Nelements+=2*nx*pmlNy;
  if (rank_y==size_y-1) Nelements+=2*nx*pmlNy;

  //pml corners
  if (rank_x==0        && rank_y==0       ) Nelements+=2*pmlNx*pmlNy;
  if (rank_x==size_x-1 && rank_y==0       ) Nelements+=2*pmlNx*pmlNy;
  if (rank_x==0        && rank_y==size_y-1) Nelements+=2*pmlNx*pmlNy;
  if (rank_x==size_x-1 && rank_y==size_y-1) Nelements+=2*pmlNx*pmlNy;


  EToV.malloc(Nelements*Nverts);
  EX.malloc(Nelements*Nverts);
  EY.malloc(Nelements*Nverts);

  elementInfo.malloc(Nelements);

  dlong e = 0;

  //interior
  for(int j=0;j<ny;++j){
    for(int i=0;i<nx;++i){

      const hlong i0 = i+offset_x + pmlNx;
      const hlong i1 = (i+1+offset_x + pmlNx)%NnX;
      const hlong j0 = j+offset_y + pmlNy;
      const hlong j1 = (j+1+offset_y + pmlNy)%NnY;

      dfloat x0 = X0 + dx*i;
      dfloat y0 = Y0 + dy*j;

      EToV[e*Nverts+0] = i0 + j0*NnX;
      EToV[e*Nverts+1] = i1 + j0*NnX;
      EToV[e*Nverts+2] = i1 + j1*NnX;

      EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;
      EX[e*Nverts+1] = x0+dx; EY[e*Nverts+1] = y0;
      EX[e*Nverts+2] = x0+dx; EY[e*Nverts+2] = y0+dy;

      elementInfo[e] = 1; // domain
      e++;

      EToV[e*Nverts+0] = i0 + j0*NnX;
      EToV[e*Nverts+1] = i1 + j1*NnX;
      EToV[e*Nverts+2] = i0 + j1*NnX;

      EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;
      EX[e*Nverts+1] = x0+dx; EY[e*Nverts+1] = y0+dy;
      EX[e*Nverts+2] = x0;    EY[e*Nverts+2] = y0+dy;

      elementInfo[e] = 1; // domain
      e++;
    }
  }

  //X- pml
  if (rank_x==0) {
    for(int j=0;j<ny;++j){
      for(int i=0;i<pmlNx;++i){

        const hlong i0 = i+offset_x;
        const hlong i1 = (i+1+offset_x)%NnX;
        const hlong j0 = j+offset_y + pmlNy;
        const hlong j1 = (j+1+offset_y + pmlNy)%NnY;

        dfloat x0 = X0-pmlWidthx + pmldx*i;
        dfloat y0 = Y0 + dy*j;

        EToV[e*Nverts+0] = i0 + j0*NnX;
        EToV[e*Nverts+1] = i1 + j0*NnX;
        EToV[e*Nverts+2] = i1 + j1*NnX;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;
        EX[e*Nverts+1] = x0+pmldx; EY[e*Nverts+1] = y0;
        EX[e*Nverts+2] = x0+pmldx; EY[e*Nverts+2] = y0+dy;

        elementInfo[e] = 100; // X pml
        e++;

        EToV[e*Nverts+0] = i0 + j0*NnX;
        EToV[e*Nverts+1] = i1 + j1*NnX;
        EToV[e*Nverts+2] = i0 + j1*NnX;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;
        EX[e*Nverts+1] = x0+pmldx; EY[e*Nverts+1] = y0+dy;
        EX[e*Nverts+2] = x0;    EY[e*Nverts+2] = y0+dy;

        elementInfo[e] = 100; // X pml
        e++;
      }
    }
  }

  //X+ pml
  if (rank_x==size_x-1) {
    for(int j=0;j<ny;++j){
      for(int i=0;i<pmlNx;++i){

        const hlong i0 = i+offset_x + pmlNx + nx;
        const hlong i1 = (i+1+offset_x + pmlNx + nx)%NnX;
        const hlong j0 = j+offset_y + pmlNy;
        const hlong j1 = (j+1+offset_y + pmlNy)%NnY;

        dfloat x0 = X0 + dimx + pmldx*i;
        dfloat y0 = Y0 + dy*j;

        EToV[e*Nverts+0] = i0 + j0*NnX;
        EToV[e*Nverts+1] = i1 + j0*NnX;
        EToV[e*Nverts+2] = i1 + j1*NnX;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;
        EX[e*Nverts+1] = x0+pmldx; EY[e*Nverts+1] = y0;
        EX[e*Nverts+2] = x0+pmldx; EY[e*Nverts+2] = y0+dy;

        elementInfo[e] = 100; // X pml
        e++;

        EToV[e*Nverts+0] = i0 + j0*NnX;
        EToV[e*Nverts+1] = i1 + j1*NnX;
        EToV[e*Nverts+2] = i0 + j1*NnX;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;
        EX[e*Nverts+1] = x0+pmldx; EY[e*Nverts+1] = y0+dy;
        EX[e*Nverts+2] = x0;    EY[e*Nverts+2] = y0+dy;

        elementInfo[e] = 100; // X pml
        e++;
      }
    }
  }

  //Y- pml
  if (rank_y==0) {
    for(int j=0;j<pmlNy;++j){
      for(int i=0;i<nx;++i){

        const hlong i0 = i+offset_x + pmlNx;
        const hlong i1 = (i+1+offset_x + pmlNx)%NnX;
        const hlong j0 = j+offset_y;
        const hlong j1 = (j+1+offset_y)%NnY;

        dfloat x0 = X0 + dx*i;
        dfloat y0 = Y0-pmlWidthy + pmldy*j;

        EToV[e*Nverts+0] = i0 + j0*NnX;
        EToV[e*Nverts+1] = i1 + j0*NnX;
        EToV[e*Nverts+2] = i1 + j1*NnX;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;
        EX[e*Nverts+1] = x0+dx; EY[e*Nverts+1] = y0;
        EX[e*Nverts+2] = x0+dx; EY[e*Nverts+2] = y0+pmldy;

        elementInfo[e] = 200; // Y pml
        e++;

        EToV[e*Nverts+0] = i0 + j0*NnX;
        EToV[e*Nverts+1] = i1 + j1*NnX;
        EToV[e*Nverts+2] = i0 + j1*NnX;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;
        EX[e*Nverts+1] = x0+dx; EY[e*Nverts+1] = y0+pmldy;
        EX[e*Nverts+2] = x0;    EY[e*Nverts+2] = y0+pmldy;

        elementInfo[e] = 200; // Y pml
        e++;
      }
    }
  }

  //Y+ pml
  if (rank_y==size_y-1) {
    for(int j=0;j<pmlNy;++j){
      for(int i=0;i<nx;++i){

        const hlong i0 = i+offset_x + pmlNx;
        const hlong i1 = (i+1+offset_x + pmlNx)%NnX;
        const hlong j0 = j+offset_y + pmlNy + ny;
        const hlong j1 = (j+1+offset_y + pmlNy + ny)%NnY;

        dfloat x0 = X0 + dx*i;
        dfloat y0 = Y0 + dimy + pmldy*j;

        EToV[e*Nverts+0] = i0 + j0*NnX;
        EToV[e*Nverts+1] = i1 + j0*NnX;
        EToV[e*Nverts+2] = i1 + j1*NnX;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;
        EX[e*Nverts+1] = x0+dx; EY[e*Nverts+1] = y0;
        EX[e*Nverts+2] = x0+dx; EY[e*Nverts+2] = y0+pmldy;

        elementInfo[e] = 200; // Y pml
        e++;

        EToV[e*Nverts+0] = i0 + j0*NnX;
        EToV[e*Nverts+1] = i1 + j1*NnX;
        EToV[e*Nverts+2] = i0 + j1*NnX;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;
        EX[e*Nverts+1] = x0+dx; EY[e*Nverts+1] = y0+pmldy;
        EX[e*Nverts+2] = x0;    EY[e*Nverts+2] = y0+pmldy;

        elementInfo[e] = 200; // Y pml
        e++;
      }
    }
  }

  //X-Y- pml
  if (rank_x==0 && rank_y==0) {
    for(int j=0;j<pmlNy;++j){
      for(int i=0;i<pmlNx;++i){

        const hlong i0 = i+offset_x;
        const hlong i1 = (i+1+offset_x)%NnX;
        const hlong j0 = j+offset_y;
        const hlong j1 = (j+1+offset_y)%NnY;

        dfloat x0 = X0-pmlWidthx + pmldx*i;
        dfloat y0 = Y0-pmlWidthy + pmldy*j;

        EToV[e*Nverts+0] = i0 + j0*NnX;
        EToV[e*Nverts+1] = i1 + j0*NnX;
        EToV[e*Nverts+2] = i1 + j1*NnX;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;
        EX[e*Nverts+1] = x0+pmldx; EY[e*Nverts+1] = y0;
        EX[e*Nverts+2] = x0+pmldx; EY[e*Nverts+2] = y0+pmldy;

        elementInfo[e] = 300; // XY pml
        e++;

        EToV[e*Nverts+0] = i0 + j0*NnX;
        EToV[e*Nverts+1] = i1 + j1*NnX;
        EToV[e*Nverts+2] = i0 + j1*NnX;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;
        EX[e*Nverts+1] = x0+pmldx; EY[e*Nverts+1] = y0+pmldy;
        EX[e*Nverts+2] = x0;    EY[e*Nverts+2] = y0+pmldy;

        elementInfo[e] = 300; // XY pml
        e++;
      }
    }
  }

  //X+Y- pml
  if (rank_x==size_x-1 && rank_y==0) {
    for(int j=0;j<pmlNy;++j){
      for(int i=0;i<pmlNx;++i){

        const hlong i0 = i+offset_x + pmlNx + nx;
        const hlong i1 = (i+1+offset_x + pmlNx + nx)%NnX;
        const hlong j0 = j+offset_y;
        const hlong j1 = (j+1+offset_y)%NnY;

        dfloat x0 = X0+dimx      + pmldx*i;
        dfloat y0 = Y0-pmlWidthy + pmldy*j;

        EToV[e*Nverts+0] = i0 + j0*NnX;
        EToV[e*Nverts+1] = i1 + j0*NnX;
        EToV[e*Nverts+2] = i1 + j1*NnX;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;
        EX[e*Nverts+1] = x0+pmldx; EY[e*Nverts+1] = y0;
        EX[e*Nverts+2] = x0+pmldx; EY[e*Nverts+2] = y0+pmldy;

        elementInfo[e] = 300; // XY pml
        e++;

        EToV[e*Nverts+0] = i0 + j0*NnX;
        EToV[e*Nverts+1] = i1 + j1*NnX;
        EToV[e*Nverts+2] = i0 + j1*NnX;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;
        EX[e*Nverts+1] = x0+pmldx; EY[e*Nverts+1] = y0+pmldy;
        EX[e*Nverts+2] = x0;    EY[e*Nverts+2] = y0+pmldy;

        elementInfo[e] = 300; // XY pml
        e++;
      }
    }
  }

  //X-Y+ pml
  if (rank_x==0 && rank_y==size_y-1) {
    for(int j=0;j<pmlNy;++j){
      for(int i=0;i<pmlNx;++i){

        const hlong i0 = i+offset_x;
        const hlong i1 = (i+1+offset_x)%NnX;
        const hlong j0 = j+offset_y + pmlNy + ny;
        const hlong j1 = (j+1+offset_y + pmlNy + ny)%NnY;

        dfloat x0 = X0-pmlWidthx + pmldx*i;
        dfloat y0 = Y0+dimy      + pmldy*j;

        EToV[e*Nverts+0] = i0 + j0*NnX;
        EToV[e*Nverts+1] = i1 + j0*NnX;
        EToV[e*Nverts+2] = i1 + j1*NnX;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;
        EX[e*Nverts+1] = x0+pmldx; EY[e*Nverts+1] = y0;
        EX[e*Nverts+2] = x0+pmldx; EY[e*Nverts+2] = y0+pmldy;

        elementInfo[e] = 300; // XY pml
        e++;

        EToV[e*Nverts+0] = i0 + j0*NnX;
        EToV[e*Nverts+1] = i1 + j1*NnX;
        EToV[e*Nverts+2] = i0 + j1*NnX;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;
        EX[e*Nverts+1] = x0+pmldx; EY[e*Nverts+1] = y0+pmldy;
        EX[e*Nverts+2] = x0;    EY[e*Nverts+2] = y0+pmldy;

        elementInfo[e] = 300; // XY pml
        e++;
      }
    }
  }

  //X+Y+ pml
  if (rank_x==size_x-1 && rank_y==size_y-1) {
    for(int j=0;j<pmlNy;++j){
      for(int i=0;i<pmlNx;++i){

        const hlong i0 = i+offset_x + pmlNx + nx;
        const hlong i1 = (i+1+offset_x + pmlNx + nx)%NnX;
        const hlong j0 = j+offset_y + pmlNy + ny;
        const hlong j1 = (j+1+offset_y + pmlNy + ny)%NnY;

        dfloat x0 = X0+dimx      + pmldx*i;
        dfloat y0 = Y0+dimy      + pmldy*j;

        EToV[e*Nverts+0] = i0 + j0*NnX;
        EToV[e*Nverts+1] = i1 + j0*NnX;
        EToV[e*Nverts+2] = i1 + j1*NnX;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;
        EX[e*Nverts+1] = x0+pmldx; EY[e*Nverts+1] = y0;
        EX[e*Nverts+2] = x0+pmldx; EY[e*Nverts+2] = y0+pmldy;

        elementInfo[e] = 300; // XY pml
        e++;

        EToV[e*Nverts+0] = i0 + j0*NnX;
        EToV[e*Nverts+1] = i1 + j1*NnX;
        EToV[e*Nverts+2] = i0 + j1*NnX;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;
        EX[e*Nverts+1] = x0+pmldx; EY[e*Nverts+1] = y0+pmldy;
        EX[e*Nverts+2] = x0;    EY[e*Nverts+2] = y0+pmldy;

        elementInfo[e] = 300; // XY pml
        e++;
      }
    }
  }

  if (boundaryFlag != -1) { //-1 reserved for periodic case
    NboundaryFaces = 2*NX + 2*NY;
    boundaryInfo.malloc(NboundaryFaces*(NfaceVertices+1));

    hlong bcnt = 0;

    //top and bottom
    for(hlong i=0;i<NX;++i){
      hlong vid1 = i +  0*NnX;
      hlong vid2 = i + NY*NnX;

      boundaryInfo[bcnt*3+0] = boundaryFlag;
      boundaryInfo[bcnt*3+1] = vid1 + 0;
      boundaryInfo[bcnt*3+2] = vid1 + 1;
      bcnt++;

      boundaryInfo[bcnt*3+0] = boundaryFlag;
      boundaryInfo[bcnt*3+1] = vid2 + 0;
      boundaryInfo[bcnt*3+2] = vid2 + 1;
      bcnt++;
    }

    //left and right
    for(hlong j=0;j<NY;++j){
      hlong vid1 =  0 + j*NnX;
      hlong vid2 = NX + j*NnX;

      boundaryInfo[bcnt*3+0] = boundaryFlag;
      boundaryInfo[bcnt*3+1] = vid1 + 0*NnX;
      boundaryInfo[bcnt*3+2] = vid1 + 1*NnX;
      bcnt++;

      boundaryInfo[bcnt*3+0] = boundaryFlag;
      boundaryInfo[bcnt*3+1] = vid2 + 0*NnX;
      boundaryInfo[bcnt*3+2] = vid2 + 1*NnX;
      bcnt++;
    }

  } else {
    NboundaryFaces = 0; // no boundaries
  }
}

} //namespace libp
