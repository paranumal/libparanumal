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

void meshTri3D::SetupPmlBox(){
  LIBP_ABORT(string("PMLBOX mesh not currently supprted for Tri3D meshes."))
}

void meshTri2D::SetupPmlBox(){

  dim = 2;
  Nverts = 3; // number of vertices per element
  Nfaces = 3;
  NfaceVertices = 2;

  // vertices on each face
  int faceVertices_[4][2] = {{0,1},{1,2},{2,0}};

  faceVertices = (int*) calloc(NfaceVertices*Nfaces, sizeof(int));
  memcpy(faceVertices, faceVertices_[0], NfaceVertices*Nfaces*sizeof(int));

  //local grid physical sizes
  dfloat DIMX, DIMY;
  settings.getSetting("BOX DIMX", DIMX);
  settings.getSetting("BOX DIMY", DIMY);

  //relative PML width
  const dfloat pmlScale = 0.2;

  //number of local elements in each dimension
  dlong nx, ny;
  settings.getSetting("BOX NX", nx);
  settings.getSetting("BOX NY", ny);

  int size_x = std::sqrt(size); //number of ranks in each dimension
  if (size_x*size_x != size)
    LIBP_ABORT(string("2D PMLBOX mesh requires a square number of ranks for now."))

  int boundaryFlag;
  settings.getSetting("BOX BOUNDARY FLAG", boundaryFlag);

  const int periodicFlag = (boundaryFlag == -1) ? 1 : 0;
  if (periodicFlag)
    LIBP_ABORT(string("Periodic boundary unsupported for PMLBOX mesh."))

  //local grid physical sizes
  dfloat dimx = DIMX/size_x;
  dfloat dimy = DIMY/size_x;

  //pml width
  dfloat pmlWidthx = pmlScale*DIMX;
  dfloat pmlWidthy = pmlScale*DIMY;

  //rank coordinates
  int rank_y = rank / size_x;
  int rank_x = rank % size_x;

  //bottom corner of physical domain
  dfloat X0 = -DIMX/2.0 + rank_x*dimx;
  dfloat Y0 = -DIMY/2.0 + rank_y*dimy;

  //global number of elements in each dimension
  hlong NX = size_x*nx;
  hlong NY = size_x*ny;

  dfloat dx = dimx/nx;
  dfloat dy = dimy/ny;

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
  if (rank_y==size_x-1) Nelements+=2*nx*pmlNy;

  //pml corners
  if (rank_x==0        && rank_y==0       ) Nelements+=2*pmlNx*pmlNy;
  if (rank_x==size_x-1 && rank_y==0       ) Nelements+=2*pmlNx*pmlNy;
  if (rank_x==0        && rank_y==size_x-1) Nelements+=2*pmlNx*pmlNy;
  if (rank_x==size_x-1 && rank_y==size_x-1) Nelements+=2*pmlNx*pmlNy;


  EToV = (hlong*) calloc(Nelements*Nverts, sizeof(hlong));
  EX = (dfloat*) calloc(Nelements*Nverts, sizeof(dfloat));
  EY = (dfloat*) calloc(Nelements*Nverts, sizeof(dfloat));

  elementInfo = (hlong*) calloc(Nelements, sizeof(hlong));

  dlong e = 0;

  //interior
  for(int j=0;j<ny;++j){
    for(int i=0;i<nx;++i){

      const hlong i0 = i+rank_x*nx + pmlNx;
      const hlong i1 = (i+1+rank_x*nx + pmlNx)%NnX;
      const hlong j0 = j+rank_y*ny + pmlNy;
      const hlong j1 = (j+1+rank_y*ny + pmlNy)%NnY;

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

        const hlong i0 = i+rank_x*nx;
        const hlong i1 = (i+1+rank_x*nx)%NnX;
        const hlong j0 = j+rank_y*ny + pmlNy;
        const hlong j1 = (j+1+rank_y*ny + pmlNy)%NnY;

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

        const hlong i0 = i+rank_x*nx + pmlNx + nx;
        const hlong i1 = (i+1+rank_x*nx + pmlNx + nx)%NnX;
        const hlong j0 = j+rank_y*ny + pmlNy;
        const hlong j1 = (j+1+rank_y*ny + pmlNy)%NnY;

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

        const hlong i0 = i+rank_x*nx + pmlNx;
        const hlong i1 = (i+1+rank_x*nx + pmlNx)%NnX;
        const hlong j0 = j+rank_y*ny;
        const hlong j1 = (j+1+rank_y*ny)%NnY;

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
  if (rank_y==size_x-1) {
    for(int j=0;j<pmlNy;++j){
      for(int i=0;i<nx;++i){

        const hlong i0 = i+rank_x*nx + pmlNx;
        const hlong i1 = (i+1+rank_x*nx + pmlNx)%NnX;
        const hlong j0 = j+rank_y*ny + pmlNy + ny;
        const hlong j1 = (j+1+rank_y*ny + pmlNy + ny)%NnY;

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

        const hlong i0 = i+rank_x*nx;
        const hlong i1 = (i+1+rank_x*nx)%NnX;
        const hlong j0 = j+rank_y*ny;
        const hlong j1 = (j+1+rank_y*ny)%NnY;

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

        const hlong i0 = i+rank_x*nx + pmlNx + nx;
        const hlong i1 = (i+1+rank_x*nx + pmlNx + nx)%NnX;
        const hlong j0 = j+rank_y*ny;
        const hlong j1 = (j+1+rank_y*ny)%NnY;

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
  if (rank_x==0 && rank_y==size_x-1) {
    for(int j=0;j<pmlNy;++j){
      for(int i=0;i<pmlNx;++i){

        const hlong i0 = i+rank_x*nx;
        const hlong i1 = (i+1+rank_x*nx)%NnX;
        const hlong j0 = j+rank_y*ny + pmlNy + ny;
        const hlong j1 = (j+1+rank_y*ny + pmlNy + ny)%NnY;

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
  if (rank_x==size_x-1 && rank_y==size_x-1) {
    for(int j=0;j<pmlNy;++j){
      for(int i=0;i<pmlNx;++i){

        const hlong i0 = i+rank_x*nx + pmlNx + nx;
        const hlong i1 = (i+1+rank_x*nx + pmlNx + nx)%NnX;
        const hlong j0 = j+rank_y*ny + pmlNy + ny;
        const hlong j1 = (j+1+rank_y*ny + pmlNy + ny)%NnY;

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
    boundaryInfo = (hlong*) calloc(NboundaryFaces*(NfaceVertices+1), sizeof(hlong));

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
    NboundaryFaces = 0;
    boundaryInfo = NULL; // no boundaries
  }
}
