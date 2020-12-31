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

void meshTri3D::SetupBox(){
  LIBP_ABORT(string("BOX mesh not currently supprted for Tri3D meshes."))
}

void meshTri2D::SetupBox(){

  dim = 2;
  Nverts = 3; // number of vertices per element
  Nfaces = 3;
  NfaceVertices = 2;

  // vertices on each face
  int faceVertices_[4][2] = {{0,1},{1,2},{2,0}};

  faceVertices = (int*) calloc(NfaceVertices*Nfaces, sizeof(int));
  memcpy(faceVertices, faceVertices_[0], NfaceVertices*Nfaces*sizeof(int));

  // find a factorization size = size_x*size_y such that
  //  size_x>=size_y and are 'close' to one another
  int size_x, size_y;
  factor2(size, size_x, size_y);

  //find our coordinates in the MPI grid such that
  // rank = rank_x + rank_y*size_x
  int rank_y = rank / size_x;
  int rank_x = rank % size_x;

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

  //grid physical sizes
  dfloat DIMX, DIMY;
  settings.getSetting("BOX DIMX", DIMX);
  settings.getSetting("BOX DIMY", DIMY);

  //element sizes
  dfloat dx = DIMX/NX;
  dfloat dy = DIMY/NY;

  dlong offset_x = rank_x*(NX/size_x) + mymin(rank_x, (NX % size_x));
  dlong offset_y = rank_y*(NY/size_y) + mymin(rank_y, (NY % size_y));

  //bottom corner of physical domain
  dfloat X0 = -DIMX/2.0 + offset_x*dx;
  dfloat Y0 = -DIMY/2.0 + offset_y*dy;

  //global number of nodes in each dimension
  hlong NnX = periodicFlag ? NX : NX+1; //lose a node when periodic (repeated node)
  hlong NnY = periodicFlag ? NY : NY+1; //lose a node when periodic (repeated node)

  // build an nx x ny x nz box grid
  Nnodes = NnX*NnY; //global node count
  Nelements = 2*nx*ny; //local

  EToV = (hlong*) calloc(Nelements*Nverts, sizeof(hlong));
  EX = (dfloat*) calloc(Nelements*Nverts, sizeof(dfloat));
  EY = (dfloat*) calloc(Nelements*Nverts, sizeof(dfloat));

  elementInfo = (hlong*) calloc(Nelements, sizeof(hlong));

  dlong e = 0;
  for(int j=0;j<ny;++j){
    for(int i=0;i<nx;++i){

      const hlong i0 = i+offset_x;
      const hlong i1 = (i+1+offset_x)%NnX;
      const hlong j0 = j+offset_y;
      const hlong j1 = (j+1+offset_y)%NnY;

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
