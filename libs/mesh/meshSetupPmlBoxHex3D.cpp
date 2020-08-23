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

static void addHex(hlong i0, hlong j0, hlong k0, hlong NnX, hlong NnY, hlong NnZ,
                    dfloat x0, dfloat y0, dfloat z0, dfloat dx, dfloat dy, dfloat dz,
                    hlong *EToV, dfloat *EX, dfloat *EY, dfloat *EZ,
                    hlong *elementInfo, int type, dlong &e);

void meshHex3D::SetupPmlBox(){

  dim = 3;
  Nverts = 8; // number of vertices per element
  Nfaces = 6;
  NfaceVertices = 4;

  // vertices on each face
  int _faceVertices[6][4] =
    {{0,1,2,3},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7},{4,5,6,7}};

  faceVertices =
    (int*) calloc(NfaceVertices*Nfaces, sizeof(int));

  memcpy(faceVertices, _faceVertices[0], NfaceVertices*Nfaces*sizeof(int));

  //local grid physical sizes
  dfloat DIMX, DIMY, DIMZ;
  settings.getSetting("BOX DIMX", DIMX);
  settings.getSetting("BOX DIMY", DIMY);
  settings.getSetting("BOX DIMZ", DIMZ);

  //relative PML width
  const dfloat pmlScale = 0.2;

  //number of local elements in each dimension
  dlong nx, ny, nz;
  settings.getSetting("BOX NX", nx);
  settings.getSetting("BOX NY", ny);
  settings.getSetting("BOX NZ", nz);

  int size_x = std::cbrt(size); //number of ranks in each dimension
  if (size_x*size_x*size_x != size)
    LIBP_ABORT(string("3D PMLBOX mesh requires a cubic number of ranks for now."))

  int boundaryFlag;
  settings.getSetting("BOX BOUNDARY FLAG", boundaryFlag);

  const int periodicFlag = (boundaryFlag == -1) ? 1 : 0;
  if (periodicFlag)
    LIBP_ABORT(string("Periodic boundary unsupported for PMLBOX mesh."))

  //local grid physical sizes
  dfloat dimx = DIMX/size_x;
  dfloat dimy = DIMY/size_x;
  dfloat dimz = DIMZ/size_x;

  //pml width
  dfloat pmlWidthx = pmlScale*DIMX;
  dfloat pmlWidthy = pmlScale*DIMY;
  dfloat pmlWidthz = pmlScale*DIMZ;

  //rank coordinates
  int rank_z = rank / (size_x*size_x);
  int rank_y = (rank - rank_z*size_x*size_x) / size_x;
  int rank_x = rank % size_x;

  //bottom corner of physical domain
  dfloat X0 = -DIMX/2.0 + rank_x*dimx;
  dfloat Y0 = -DIMY/2.0 + rank_y*dimy;
  dfloat Z0 = -DIMZ/2.0 + rank_z*dimz;

  //global number of elements in each dimension
  hlong NX = size_x*nx;
  hlong NY = size_x*ny;
  hlong NZ = size_x*nz;

  dfloat dx = dimx/nx;
  dfloat dy = dimy/ny;
  dfloat dz = dimz/nz;

  //try and have roughly the same resolution in the pml region
  dlong pmlNx, pmlNy, pmlNz;
  pmlNx = ceil(pmlWidthx/dx);
  pmlNy = ceil(pmlWidthy/dy);
  pmlNz = ceil(pmlWidthz/dz);

  dfloat pmldx, pmldy, pmldz;
  pmldx = pmlWidthx/pmlNx;
  pmldy = pmlWidthy/pmlNy;
  pmldz = pmlWidthz/pmlNz;

  NX += 2*pmlNx;
  NY += 2*pmlNy;
  NZ += 2*pmlNz;

  //global number of nodes in each dimension
  hlong NnX = periodicFlag ? NX : NX+1; //lose a node when periodic (repeated node)
  hlong NnY = periodicFlag ? NY : NY+1; //lose a node when periodic (repeated node)
  hlong NnZ = periodicFlag ? NZ : NZ+1; //lose a node when periodic (repeated node)

  // build an nx x ny x nz box grid
  Nnodes = NnX*NnY*NnZ; //global node count
  Nelements = nx*ny*nz; //local

  //pml faces
  if (rank_x==0)        Nelements+=ny*nz*pmlNx;
  if (rank_x==size_x-1) Nelements+=ny*nz*pmlNx;
  if (rank_y==0)        Nelements+=nx*nz*pmlNy;
  if (rank_y==size_x-1) Nelements+=nx*nz*pmlNy;
  if (rank_z==0)        Nelements+=nx*ny*pmlNz;
  if (rank_z==size_x-1) Nelements+=nx*ny*pmlNz;

  // //pml edges
  if (rank_x==0        && rank_y==0       ) Nelements+=pmlNx*pmlNy*nz;
  if (rank_x==size_x-1 && rank_y==0       ) Nelements+=pmlNx*pmlNy*nz;
  if (rank_x==0        && rank_y==size_x-1) Nelements+=pmlNx*pmlNy*nz;
  if (rank_x==size_x-1 && rank_y==size_x-1) Nelements+=pmlNx*pmlNy*nz;
  if (rank_x==0        && rank_z==0       ) Nelements+=pmlNx*pmlNz*ny;
  if (rank_x==size_x-1 && rank_z==0       ) Nelements+=pmlNx*pmlNz*ny;
  if (rank_x==0        && rank_z==size_x-1) Nelements+=pmlNx*pmlNz*ny;
  if (rank_x==size_x-1 && rank_z==size_x-1) Nelements+=pmlNx*pmlNz*ny;
  if (rank_y==0        && rank_z==0       ) Nelements+=pmlNy*pmlNz*nx;
  if (rank_y==size_x-1 && rank_z==0       ) Nelements+=pmlNy*pmlNz*nx;
  if (rank_y==0        && rank_z==size_x-1) Nelements+=pmlNy*pmlNz*nx;
  if (rank_y==size_x-1 && rank_z==size_x-1) Nelements+=pmlNy*pmlNz*nx;

  //pml corners
  if (rank_x==0        && rank_y==0        && rank_z==0       ) Nelements+=pmlNx*pmlNy*pmlNz;
  if (rank_x==size_x-1 && rank_y==0        && rank_z==0       ) Nelements+=pmlNx*pmlNy*pmlNz;
  if (rank_x==0        && rank_y==size_x-1 && rank_z==0       ) Nelements+=pmlNx*pmlNy*pmlNz;
  if (rank_x==size_x-1 && rank_y==size_x-1 && rank_z==0       ) Nelements+=pmlNx*pmlNy*pmlNz;
  if (rank_x==0        && rank_y==0        && rank_z==size_x-1) Nelements+=pmlNx*pmlNy*pmlNz;
  if (rank_x==size_x-1 && rank_y==0        && rank_z==size_x-1) Nelements+=pmlNx*pmlNy*pmlNz;
  if (rank_x==0        && rank_y==size_x-1 && rank_z==size_x-1) Nelements+=pmlNx*pmlNy*pmlNz;
  if (rank_x==size_x-1 && rank_y==size_x-1 && rank_z==size_x-1) Nelements+=pmlNx*pmlNy*pmlNz;


  EToV = (hlong*) calloc(Nelements*Nverts, sizeof(hlong));
  EX = (dfloat*) calloc(Nelements*Nverts, sizeof(dfloat));
  EY = (dfloat*) calloc(Nelements*Nverts, sizeof(dfloat));
  EZ = (dfloat*) calloc(Nelements*Nverts, sizeof(dfloat));

  elementInfo = (hlong*) calloc(Nelements, sizeof(hlong));

  dlong e = 0;
  for(int k=0;k<nz;++k){
    for(int j=0;j<ny;++j){
      for(int i=0;i<nx;++i){

        const hlong i0 = i+rank_x*nx + pmlNx;
        const hlong j0 = j+rank_y*ny + pmlNy;
        const hlong k0 = k+rank_z*nz + pmlNz;

        dfloat x0 = X0 + dx*i;
        dfloat y0 = Y0 + dy*j;
        dfloat z0 = Z0 + dz*k;

        int type = 1; // domain

        addHex(i0, j0, k0, NnX, NnY, NnZ,
                x0, y0, z0, dx, dy, dz,
                EToV, EX, EY, EZ,
                elementInfo, type, e);
      }
    }
  }

  //X- pml
  if (rank_x==0) {
    for(int k=0;k<nz;++k){
      for(int j=0;j<ny;++j){
        for(int i=0;i<pmlNx;++i){

          const hlong i0 = i+rank_x*nx;
          const hlong j0 = j+rank_y*ny + pmlNy;
          const hlong k0 = k+rank_z*nz + pmlNz;

          dfloat x0 = X0-pmlWidthx + pmldx*i;
          dfloat y0 = Y0 + dy*j;
          dfloat z0 = Z0 + dz*k;

          int type = 100; // X pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, pmldx, dy, dz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }

  //X+ pml
  if (rank_x==size_x-1) {
    for(int k=0;k<nz;++k){
      for(int j=0;j<ny;++j){
        for(int i=0;i<pmlNx;++i){

          const hlong i0 = i+rank_x*nx + pmlNx + nx;
          const hlong j0 = j+rank_y*ny + pmlNy;
          const hlong k0 = k+rank_z*nz + pmlNz;

          dfloat x0 = X0 + dimx + pmldx*i;
          dfloat y0 = Y0 + dy*j;
          dfloat z0 = Z0 + dz*k;

          int type = 100; // X pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, pmldx, dy, dz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }

  //Y- pml
  if (rank_y==0) {
    for(int k=0;k<nz;++k){
      for(int j=0;j<pmlNy;++j){
        for(int i=0;i<nx;++i){

          const hlong i0 = i+rank_x*nx + pmlNx;
          const hlong j0 = j+rank_y*ny;
          const hlong k0 = k+rank_z*nz + pmlNz;

          dfloat x0 = X0 + dx*i;
          dfloat y0 = Y0-pmlWidthy + pmldy*j;
          dfloat z0 = Z0 + dz*k;

          int type = 200; // Y pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, dx, pmldy, dz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }

  //Y+ pml
  if (rank_y==size_x-1) {
    for(int k=0;k<nz;++k){
      for(int j=0;j<pmlNy;++j){
        for(int i=0;i<nx;++i){

          const hlong i0 = i+rank_x*nx + pmlNx;
          const hlong j0 = j+rank_y*ny + pmlNy + ny;
          const hlong k0 = k+rank_z*nz + pmlNz;

          dfloat x0 = X0 + dx*i;
          dfloat y0 = Y0 + dimy + pmldy*j;
          dfloat z0 = Z0 + dz*k;

          int type = 200; // Y pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, dx, pmldy, dz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }

  //Z- pml
  if (rank_z==0) {
    for(int k=0;k<pmlNz;++k){
      for(int j=0;j<ny;++j){
        for(int i=0;i<nx;++i){

          const hlong i0 = i+rank_x*nx + pmlNx;
          const hlong j0 = j+rank_y*ny + pmlNy;
          const hlong k0 = k+rank_z*nz;

          dfloat x0 = X0 + dx*i;
          dfloat y0 = Y0 + dy*j;
          dfloat z0 = Z0-pmlWidthz + pmldz*k;

          int type = 400; // Z pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, dx, dy, pmldz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }

  //Z+ pml
  if (rank_z==size_x-1) {
    for(int k=0;k<pmlNz;++k){
      for(int j=0;j<ny;++j){
        for(int i=0;i<nx;++i){

          const hlong i0 = i+rank_x*nx + pmlNx;
          const hlong j0 = j+rank_y*ny + pmlNy;
          const hlong k0 = k+rank_z*nz + pmlNz + nz;

          dfloat x0 = X0 + dx*i;
          dfloat y0 = Y0 + dy*j;
          dfloat z0 = Z0 + dimz + pmldz*k;

          int type = 400; // Z pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, dx, dy, pmldz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }

  //X-Y- pml
  if (rank_x==0 && rank_y==0) {
    for(int k=0;k<nz;++k){
      for(int j=0;j<pmlNy;++j){
        for(int i=0;i<pmlNx;++i){

          const hlong i0 = i+rank_x*nx;
          const hlong j0 = j+rank_y*ny;
          const hlong k0 = k+rank_z*nz + pmlNz;

          dfloat x0 = X0-pmlWidthx + pmldx*i;
          dfloat y0 = Y0-pmlWidthy + pmldy*j;
          dfloat z0 = Z0 + dz*k;

          int type = 300; // XY pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, pmldx, pmldy, dz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }
  //X+Y- pml
  if (rank_x==size_x-1 && rank_y==0) {
    for(int k=0;k<nz;++k){
      for(int j=0;j<pmlNy;++j){
        for(int i=0;i<pmlNx;++i){

          const hlong i0 = i+rank_x*nx + pmlNx + nx;
          const hlong j0 = j+rank_y*ny;
          const hlong k0 = k+rank_z*nz + pmlNz;

          dfloat x0 = X0+dimx + pmldx*i;
          dfloat y0 = Y0-pmlWidthy + pmldy*j;
          dfloat z0 = Z0 + dz*k;

          int type = 300; // XY pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, pmldx, pmldy, dz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }
  //X-Y+ pml
  if (rank_x==0 && rank_y==size_x-1) {
    for(int k=0;k<nz;++k){
      for(int j=0;j<pmlNy;++j){
        for(int i=0;i<pmlNx;++i){

          const hlong i0 = i+rank_x*nx;
          const hlong j0 = j+rank_y*ny + pmlNy + ny;
          const hlong k0 = k+rank_z*nz + pmlNz;

          dfloat x0 = X0-pmlWidthx + pmldx*i;
          dfloat y0 = Y0+dimy + pmldy*j;
          dfloat z0 = Z0 + dz*k;

          int type = 300; // XY pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, pmldx, pmldy, dz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }
  //X+Y+ pml
  if (rank_x==size_x-1 && rank_y==size_x-1) {
    for(int k=0;k<nz;++k){
      for(int j=0;j<pmlNy;++j){
        for(int i=0;i<pmlNx;++i){

          const hlong i0 = i+rank_x*nx + pmlNx + nx;
          const hlong j0 = j+rank_y*ny + pmlNy + ny;
          const hlong k0 = k+rank_z*nz + pmlNz;

          dfloat x0 = X0+dimx + pmldx*i;
          dfloat y0 = Y0+dimy + pmldy*j;
          dfloat z0 = Z0 + dz*k;

          int type = 300; // XY pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, pmldx, pmldy, dz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }

  //X-Z- pml
  if (rank_x==0 && rank_z==0) {
    for(int k=0;k<pmlNz;++k){
      for(int j=0;j<ny;++j){
        for(int i=0;i<pmlNx;++i){

          const hlong i0 = i+rank_x*nx;
          const hlong j0 = j+rank_y*ny + pmlNy;
          const hlong k0 = k+rank_z*nz;

          dfloat x0 = X0-pmlWidthx + pmldx*i;
          dfloat y0 = Y0 + dy*j;
          dfloat z0 = Z0-pmlWidthz + pmldz*k;

          int type = 500; // XZ pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, pmldx, dy, pmldz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }
  //X+Z- pml
  if (rank_x==size_x-1 && rank_z==0) {
    for(int k=0;k<pmlNz;++k){
      for(int j=0;j<ny;++j){
        for(int i=0;i<pmlNx;++i){

          const hlong i0 = i+rank_x*nx + pmlNx + nx;
          const hlong j0 = j+rank_y*ny + pmlNy;
          const hlong k0 = k+rank_z*nz;

          dfloat x0 = X0+dimx + pmldx*i;
          dfloat y0 = Y0 + dy*j;
          dfloat z0 = Z0-pmlWidthz + pmldz*k;

          int type = 500; // XZ pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, pmldx, dy, pmldz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }
  //X-Z+ pml
  if (rank_x==0 && rank_z==size_x-1) {
    for(int k=0;k<pmlNz;++k){
      for(int j=0;j<ny;++j){
        for(int i=0;i<pmlNx;++i){

          const hlong i0 = i+rank_x*nx;
          const hlong j0 = j+rank_y*ny + pmlNy;
          const hlong k0 = k+rank_z*nz + pmlNz + nz;

          dfloat x0 = X0-pmlWidthx + pmldx*i;
          dfloat y0 = Y0 + dy*j;
          dfloat z0 = Z0+dimz + pmldz*k;

          int type = 500; // XZ pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, pmldx, dy, pmldz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }
  //X+Z+ pml
  if (rank_x==size_x-1 && rank_z==size_x-1) {
    for(int k=0;k<pmlNz;++k){
      for(int j=0;j<ny;++j){
        for(int i=0;i<pmlNx;++i){

          const hlong i0 = i+rank_x*nx + pmlNx + nx;
          const hlong j0 = j+rank_y*ny + pmlNy;
          const hlong k0 = k+rank_z*nz + pmlNz + nz;

          dfloat x0 = X0+dimx + pmldx*i;
          dfloat y0 = Y0 + dy*j;
          dfloat z0 = Z0+dimz + pmldz*k;

          int type = 500; // XZ pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, pmldx, dy, pmldz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }

  //Y-Z- pml
  if (rank_x==0 && rank_z==0) {
    for(int k=0;k<pmlNz;++k){
      for(int j=0;j<pmlNy;++j){
        for(int i=0;i<nx;++i){

          const hlong i0 = i+rank_x*nx + pmlNx;
          const hlong j0 = j+rank_y*ny;
          const hlong k0 = k+rank_z*nz;

          dfloat x0 = X0 + dx*i;
          dfloat y0 = Y0-pmlWidthy + pmldy*j;
          dfloat z0 = Z0-pmlWidthz + pmldz*k;

          int type = 600; // YZ pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, dx, pmldy, pmldz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }
  //Y+Z- pml
  if (rank_y==size_x-1 && rank_z==0) {
    for(int k=0;k<pmlNz;++k){
      for(int j=0;j<pmlNy;++j){
        for(int i=0;i<nx;++i){

          const hlong i0 = i+rank_x*nx + pmlNx;
          const hlong j0 = j+rank_y*ny + pmlNy + ny;
          const hlong k0 = k+rank_z*nz;

          dfloat x0 = X0 + dx*i;
          dfloat y0 = Y0+dimy + pmldy*j;
          dfloat z0 = Z0-pmlWidthz + pmldz*k;

          int type = 600; // YZ pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, dx, pmldy, pmldz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }
  //Y-Z+ pml
  if (rank_y==0 && rank_z==size_x-1) {
    for(int k=0;k<pmlNz;++k){
      for(int j=0;j<pmlNy;++j){
        for(int i=0;i<nx;++i){

          const hlong i0 = i+rank_x*nx + pmlNx;
          const hlong j0 = j+rank_y*ny;
          const hlong k0 = k+rank_z*nz + pmlNz + nz;

          dfloat x0 = X0 + dx*i;
          dfloat y0 = Y0-pmlWidthy + pmldy*j;
          dfloat z0 = Z0+dimz + pmldz*k;

          int type = 600; // YZ pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, dx, pmldy, pmldz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }
  //Y+Z+ pml
  if (rank_y==size_x-1 && rank_z==size_x-1) {
    for(int k=0;k<pmlNz;++k){
      for(int j=0;j<pmlNy;++j){
        for(int i=0;i<nx;++i){

          const hlong i0 = i+rank_x*nx + pmlNx;
          const hlong j0 = j+rank_y*ny + pmlNy + ny;
          const hlong k0 = k+rank_z*nz + pmlNz + nz;

          dfloat x0 = X0 + dx*i;
          dfloat y0 = Y0+dimy + pmldy*j;
          dfloat z0 = Z0+dimz + pmldz*k;

          int type = 600; // YZ pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, dx, pmldy, pmldz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }

  //X-Y-Z- pml
  if (rank_x==0 && rank_y==0 && rank_z==0) {
    for(int k=0;k<pmlNz;++k){
      for(int j=0;j<pmlNy;++j){
        for(int i=0;i<pmlNx;++i){

          const hlong i0 = i+rank_x*nx;
          const hlong j0 = j+rank_y*ny;
          const hlong k0 = k+rank_z*nz;

          dfloat x0 = X0-pmlWidthx + pmldx*i;
          dfloat y0 = Y0-pmlWidthy + pmldy*j;
          dfloat z0 = Z0-pmlWidthz + pmldz*k;

          int type = 700; // XYZ pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, pmldx, pmldy, pmldz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }
  //X+Y-Z- pml
  if (rank_x==size_x-1 && rank_y==0 && rank_z==0) {
    for(int k=0;k<pmlNz;++k){
      for(int j=0;j<pmlNy;++j){
        for(int i=0;i<pmlNx;++i){

          const hlong i0 = i+rank_x*nx + pmlNx + nx;
          const hlong j0 = j+rank_y*ny;
          const hlong k0 = k+rank_z*nz;

          dfloat x0 = X0+dimx + pmldx*i;
          dfloat y0 = Y0-pmlWidthy + pmldy*j;
          dfloat z0 = Z0-pmlWidthz + pmldz*k;

          int type = 700; // XYZ pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, pmldx, pmldy, pmldz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }
  //X-Y+Z- pml
  if (rank_x==0 && rank_y==size_x-1 && rank_z==0) {
    for(int k=0;k<pmlNz;++k){
      for(int j=0;j<pmlNy;++j){
        for(int i=0;i<pmlNx;++i){

          const hlong i0 = i+rank_x*nx;
          const hlong j0 = j+rank_y*ny + pmlNy + ny;
          const hlong k0 = k+rank_z*nz;

          dfloat x0 = X0-pmlWidthx + pmldx*i;
          dfloat y0 = Y0+dimy + pmldy*j;
          dfloat z0 = Z0-pmlWidthz + pmldz*k;

          int type = 700; // XYZ pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, pmldx, pmldy, pmldz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }
  //X+Y+Z- pml
  if (rank_x==size_x-1 && rank_y==size_x-1 && rank_z==0) {
    for(int k=0;k<pmlNz;++k){
      for(int j=0;j<pmlNy;++j){
        for(int i=0;i<pmlNx;++i){

          const hlong i0 = i+rank_x*nx + pmlNx + nx;
          const hlong j0 = j+rank_y*ny + pmlNy + ny;
          const hlong k0 = k+rank_z*nz;

          dfloat x0 = X0+dimx + pmldx*i;
          dfloat y0 = Y0+dimy + pmldy*j;
          dfloat z0 = Z0-pmlWidthz + pmldz*k;

          int type = 700; // XYZ pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, pmldx, pmldy, pmldz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }
  //X-Y-Z+ pml
  if (rank_x==0 && rank_y==0 && rank_z==size_x-1) {
    for(int k=0;k<pmlNz;++k){
      for(int j=0;j<pmlNy;++j){
        for(int i=0;i<pmlNx;++i){

          const hlong i0 = i+rank_x*nx;
          const hlong j0 = j+rank_y*ny;
          const hlong k0 = k+rank_z*nz + pmlNz + nz;

          dfloat x0 = X0-pmlWidthx + pmldx*i;
          dfloat y0 = Y0-pmlWidthy + pmldy*j;
          dfloat z0 = Z0+dimz + pmldz*k;

          int type = 700; // XYZ pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, pmldx, pmldy, pmldz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }
  //X+Y-Z+ pml
  if (rank_x==size_x-1 && rank_y==0 && rank_z==size_x-1) {
    for(int k=0;k<pmlNz;++k){
      for(int j=0;j<pmlNy;++j){
        for(int i=0;i<pmlNx;++i){

          const hlong i0 = i+rank_x*nx + pmlNx + nx;
          const hlong j0 = j+rank_y*ny;
          const hlong k0 = k+rank_z*nz + pmlNz + nz;

          dfloat x0 = X0+dimx + pmldx*i;
          dfloat y0 = Y0-pmlWidthy + pmldy*j;
          dfloat z0 = Z0+dimz + pmldz*k;

          int type = 700; // XYZ pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, pmldx, pmldy, pmldz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }
  //X-Y+Z+ pml
  if (rank_x==0 && rank_y==size_x-1 && rank_z==size_x-1) {
    for(int k=0;k<pmlNz;++k){
      for(int j=0;j<pmlNy;++j){
        for(int i=0;i<pmlNx;++i){

          const hlong i0 = i+rank_x*nx;
          const hlong j0 = j+rank_y*ny + pmlNy + ny;
          const hlong k0 = k+rank_z*nz + pmlNz + nz;

          dfloat x0 = X0-pmlWidthx + pmldx*i;
          dfloat y0 = Y0+dimy + pmldy*j;
          dfloat z0 = Z0+dimz + pmldz*k;

          int type = 700; // XYZ pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, pmldx, pmldy, pmldz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }
  //X+Y+Z+ pml
  if (rank_x==size_x-1 && rank_y==size_x-1 && rank_z==size_x-1) {
    for(int k=0;k<pmlNz;++k){
      for(int j=0;j<pmlNy;++j){
        for(int i=0;i<pmlNx;++i){

          const hlong i0 = i+rank_x*nx + pmlNx + nx;
          const hlong j0 = j+rank_y*ny + pmlNy + ny;
          const hlong k0 = k+rank_z*nz + pmlNz + nz;

          dfloat x0 = X0+dimx + pmldx*i;
          dfloat y0 = Y0+dimy + pmldy*j;
          dfloat z0 = Z0+dimz + pmldz*k;

          int type = 700; // XYZ pml

          addHex(i0, j0, k0, NnX, NnY, NnZ,
                  x0, y0, z0, pmldx, pmldy, pmldz,
                  EToV, EX, EY, EZ,
                  elementInfo, type, e);
        }
      }
    }
  }



  if (boundaryFlag != -1) { //-1 reserved for periodic case
    NboundaryFaces = 2*NX*NY + 2*NX*NZ + 2*NY*NZ;
    boundaryInfo = (hlong*) calloc(NboundaryFaces*(NfaceVertices+1), sizeof(hlong));

    hlong bcnt = 0;

    //top and bottom
    for(hlong j=0;j<NY;++j){
      for(hlong i=0;i<NX;++i){
        hlong vid1 = i + j*NnX +  0*NnX*NnY;
        hlong vid2 = i + j*NnX + NZ*NnX*NnY;

        boundaryInfo[bcnt*5+0] = boundaryFlag;
        boundaryInfo[bcnt*5+1] = vid1 + 0;
        boundaryInfo[bcnt*5+2] = vid1 + 1;
        boundaryInfo[bcnt*5+3] = vid1 + 1 + NnX;
        boundaryInfo[bcnt*5+4] = vid1 + 0 + NnX;
        bcnt++;

        boundaryInfo[bcnt*5+0] = boundaryFlag;
        boundaryInfo[bcnt*5+1] = vid2 + 0;
        boundaryInfo[bcnt*5+2] = vid2 + 1;
        boundaryInfo[bcnt*5+3] = vid2 + 1 + NnX;
        boundaryInfo[bcnt*5+4] = vid2 + 0 + NnX;
        bcnt++;
      }
    }
    //front and back
    for(hlong k=0;k<NZ;++k){
      for(hlong i=0;i<NX;++i){
        hlong vid1 = i +  0*NnX + k*NnX*NnY;
        hlong vid2 = i + NY*NnX + k*NnX*NnY;

        boundaryInfo[bcnt*5+0] = boundaryFlag;
        boundaryInfo[bcnt*5+1] = vid1 + 0;
        boundaryInfo[bcnt*5+2] = vid1 + 1;
        boundaryInfo[bcnt*5+3] = vid1 + 1 + NnX*NnY;
        boundaryInfo[bcnt*5+4] = vid1 + 0 + NnX*NnY;
        bcnt++;

        boundaryInfo[bcnt*5+0] = boundaryFlag;
        boundaryInfo[bcnt*5+1] = vid2 + 0;
        boundaryInfo[bcnt*5+2] = vid2 + 1;
        boundaryInfo[bcnt*5+3] = vid2 + 1 + NnX*NnY;
        boundaryInfo[bcnt*5+4] = vid2 + 0 + NnX*NnY;
        bcnt++;
      }
    }
    //left and right
    for(hlong k=0;k<NZ;++k){
      for(hlong j=0;j<NY;++j){
        hlong vid1 =  0 + j*NnX + k*NnX*NnY;
        hlong vid2 = NX + j*NnX + k*NnX*NnY;

        boundaryInfo[bcnt*5+0] = boundaryFlag;
        boundaryInfo[bcnt*5+1] = vid1 + 0*NnX;
        boundaryInfo[bcnt*5+2] = vid1 + 1*NnX;
        boundaryInfo[bcnt*5+3] = vid1 + 1*NnX + NnX*NnY;
        boundaryInfo[bcnt*5+4] = vid1 + 0*NnX + NnX*NnY;
        bcnt++;

        boundaryInfo[bcnt*5+0] = boundaryFlag;
        boundaryInfo[bcnt*5+1] = vid2 + 0*NnX;
        boundaryInfo[bcnt*5+2] = vid2 + 1*NnX;
        boundaryInfo[bcnt*5+3] = vid2 + 1*NnX + NnX*NnY;
        boundaryInfo[bcnt*5+4] = vid2 + 0*NnX + NnX*NnY;
        bcnt++;
      }
    }

  } else {
    NboundaryFaces = 0;
    boundaryInfo = NULL; // no boundaries
  }
}

static void addHex(hlong i0, hlong j0, hlong k0, hlong NnX, hlong NnY, hlong NnZ,
                    dfloat x0, dfloat y0, dfloat z0, dfloat dx, dfloat dy, dfloat dz,
                    hlong *EToV, dfloat *EX, dfloat *EY, dfloat *EZ,
                    hlong *elementInfo, int type, dlong &e) {

  const hlong i1 = (i0+1)%NnX;
  const hlong j1 = (j0+1)%NnY;
  const hlong k1 = (k0+1)%NnZ;

  const int Nverts = 8;

  EToV[e*Nverts+0] = i0 + j0*NnX + k0*NnX*NnY;
  EToV[e*Nverts+1] = i1 + j0*NnX + k0*NnX*NnY;
  EToV[e*Nverts+2] = i1 + j1*NnX + k0*NnX*NnY;
  EToV[e*Nverts+3] = i0 + j1*NnX + k0*NnX*NnY;

  EToV[e*Nverts+4] = i0 + j0*NnX + k1*NnX*NnY;
  EToV[e*Nverts+5] = i1 + j0*NnX + k1*NnX*NnY;
  EToV[e*Nverts+6] = i1 + j1*NnX + k1*NnX*NnY;
  EToV[e*Nverts+7] = i0 + j1*NnX + k1*NnX*NnY;

  dfloat *ex = EX+e*Nverts;
  dfloat *ey = EY+e*Nverts;
  dfloat *ez = EZ+e*Nverts;

  ex[0] = x0;    ey[0] = y0;    ez[0] = z0;
  ex[1] = x0+dx; ey[1] = y0;    ez[1] = z0;
  ex[2] = x0+dx; ey[2] = y0+dy; ez[2] = z0;
  ex[3] = x0;    ey[3] = y0+dy; ez[3] = z0;

  ex[4] = x0;    ey[4] = y0;    ez[4] = z0+dz;
  ex[5] = x0+dx; ey[5] = y0;    ez[5] = z0+dz;
  ex[6] = x0+dx; ey[6] = y0+dy; ez[6] = z0+dz;
  ex[7] = x0;    ey[7] = y0+dy; ez[7] = z0+dz;

  elementInfo[e] = type;
  e++;
}