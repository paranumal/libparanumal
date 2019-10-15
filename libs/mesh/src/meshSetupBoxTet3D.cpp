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

void meshTet3D::SetupBox(){

  dim = 3;
  Nverts = 4; // number of vertices per element
  Nfaces = 4;
  NfaceVertices = 3;

  // vertices on each face
  int faceVertices_[4][3] = {{0,1,2},{0,1,3},{1,2,3},{2,0,3}};

  faceVertices = (int*) calloc(NfaceVertices*Nfaces, sizeof(int));
  memcpy(faceVertices, faceVertices_[0], 12*sizeof(int));

  //local grid physical sizes
  dfloat DIMX, DIMY, DIMZ;
  settings.getSetting("BOX DIMX", DIMX);
  settings.getSetting("BOX DIMY", DIMY);
  settings.getSetting("BOX DIMZ", DIMZ);

  //number of local elements in each dimension
  dlong nx, ny, nz;
  settings.getSetting("BOX NX", nx);
  settings.getSetting("BOX NY", ny);
  settings.getSetting("BOX NZ", nz);

  int size_x = std::cbrt(size); //number of ranks in each dimension
  if (size_x*size_x*size_x != size)
    LIBP_ABORT(string("3D BOX mesh requires a cubic number of ranks for now."))

  int boundaryFlag;
  settings.getSetting("BOX BOUNDARY FLAG", boundaryFlag);

  const int periodicFlag = (boundaryFlag == -1) ? 1 : 0;

  //local grid physical sizes
  dfloat dimx = DIMX/size_x;
  dfloat dimy = DIMY/size_x;
  dfloat dimz = DIMZ/size_x;

  //rank coordinates
  int rank_z = rank / (size_x*size_x);
  int rank_y = (rank - rank_z*size_x*size_x) / size_x;
  int rank_x = rank % size_x;

  // //bottom corner of physical domain
  // dfloat X0 = -DIMX/2.0 + rank_x*dimx;
  // dfloat Y0 = -DIMY/2.0 + rank_y*dimy;
  // dfloat Z0 = -DIMZ/2.0 + rank_z*dimz;

  dfloat domx = 0.0; // -DIMX/2.0 ; 
  dfloat domy = 0.0; // -DIMY/2.0 ; 
  dfloat domz = 0.0; // -DIMZ/2.0 ; 

  settings.getSetting("BOX DOMX", domx);
  settings.getSetting("BOX DOMY", domy);
  settings.getSetting("BOX DOMZ", domz);

  //bottom corner of physical domain
  dfloat X0 = domx + rank_x*dimx;
  dfloat Y0 = domy + rank_y*dimy;
  dfloat Z0 = domz + rank_z*dimz;

  //global number of elements in each dimension
  hlong NX = size_x*nx;
  hlong NY = size_x*ny;
  hlong NZ = size_x*nz;

  //global number of nodes in each dimension
  hlong NnX = periodicFlag ? NX : NX+1; //lose a node when periodic (repeated node)
  hlong NnY = periodicFlag ? NY : NY+1; //lose a node when periodic (repeated node)
  hlong NnZ = periodicFlag ? NZ : NZ+1; //lose a node when periodic (repeated node)

  // build an nx x ny x nz box grid
  Nnodes = NnX*NnY*NnZ; //global node count
  Nelements = 6*nx*ny*nz; //local element count (each cube divided into 6 tets)

  EToV = (hlong*) calloc(Nelements*Nverts, sizeof(hlong));
  EX = (dfloat*) calloc(Nelements*Nverts, sizeof(dfloat));
  EY = (dfloat*) calloc(Nelements*Nverts, sizeof(dfloat));
  EZ = (dfloat*) calloc(Nelements*Nverts, sizeof(dfloat));

  elementInfo = (hlong*) calloc(Nelements, sizeof(hlong));

  dlong e = 0;
  dfloat dx = dimx/nx;
  dfloat dy = dimy/ny;
  dfloat dz = dimz/nz;
  for(int k=0;k<nz;++k){
    for(int j=0;j<ny;++j){
      for(int i=0;i<nx;++i){

        const hlong i0 = i+rank_x*nx;
        const hlong i1 = (i+1+rank_x*nx)%NnX;
        const hlong j0 = j+rank_y*ny;
        const hlong j1 = (j+1+rank_y*ny)%NnY;
        const hlong k0 = k+rank_z*nz;
        const hlong k1 = (k+1+rank_z*nz)%NnZ;

        dfloat x0 = X0 + dx*i;
        dfloat y0 = Y0 + dy*j;
        dfloat z0 = Z0 + dz*k;

        //tet 1 (0,3,2,7)
        EToV[e*Nverts+0] = i0 + j0*NnX + k0*NnX*NnY;
        EToV[e*Nverts+1] = i1 + j1*NnX + k0*NnX*NnY;
        EToV[e*Nverts+2] = i0 + j1*NnX + k0*NnX*NnY;
        EToV[e*Nverts+3] = i1 + j1*NnX + k1*NnX*NnY;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;    EZ[e*Nverts+0] = z0;
        EX[e*Nverts+1] = x0+dx; EY[e*Nverts+1] = y0+dy; EZ[e*Nverts+1] = z0;
        EX[e*Nverts+2] = x0;    EY[e*Nverts+2] = y0+dy; EZ[e*Nverts+2] = z0;
        EX[e*Nverts+3] = x0+dx; EY[e*Nverts+3] = y0+dy; EZ[e*Nverts+3] = z0+dz;

        elementInfo[e] = 1; // domain
        e++;

        //tet 2 (0,1,3,7)
        EToV[e*Nverts+0] = i0 + j0*NnX + k0*NnX*NnY;
        EToV[e*Nverts+1] = i1 + j0*NnX + k0*NnX*NnY;
        EToV[e*Nverts+2] = i1 + j1*NnX + k0*NnX*NnY;
        EToV[e*Nverts+3] = i1 + j1*NnX + k1*NnX*NnY;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;    EZ[e*Nverts+0] = z0;
        EX[e*Nverts+1] = x0+dx; EY[e*Nverts+1] = y0;    EZ[e*Nverts+1] = z0;
        EX[e*Nverts+2] = x0+dx; EY[e*Nverts+2] = y0+dy; EZ[e*Nverts+2] = z0;
        EX[e*Nverts+3] = x0+dx; EY[e*Nverts+3] = y0+dy; EZ[e*Nverts+3] = z0+dz;

        elementInfo[e] = 1; // domain
        e++;

        //tet 3 (0,2,6,7)
        EToV[e*Nverts+0] = i0 + j0*NnX + k0*NnX*NnY;
        EToV[e*Nverts+1] = i0 + j1*NnX + k0*NnX*NnY;
        EToV[e*Nverts+2] = i0 + j1*NnX + k1*NnX*NnY;
        EToV[e*Nverts+3] = i1 + j1*NnX + k1*NnX*NnY;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;    EZ[e*Nverts+0] = z0;
        EX[e*Nverts+1] = x0;    EY[e*Nverts+1] = y0+dy; EZ[e*Nverts+1] = z0;
        EX[e*Nverts+2] = x0;    EY[e*Nverts+2] = y0+dy; EZ[e*Nverts+2] = z0+dz;
        EX[e*Nverts+3] = x0+dx; EY[e*Nverts+3] = y0+dy; EZ[e*Nverts+3] = z0+dz;

        elementInfo[e] = 1; // domain
        e++;

        //tet 4 (0,6,4,7)
        EToV[e*Nverts+0] = i0 + j0*NnX + k0*NnX*NnY;
        EToV[e*Nverts+1] = i0 + j1*NnX + k1*NnX*NnY;
        EToV[e*Nverts+2] = i0 + j0*NnX + k1*NnX*NnY;
        EToV[e*Nverts+3] = i1 + j1*NnX + k1*NnX*NnY;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;    EZ[e*Nverts+0] = z0;
        EX[e*Nverts+1] = x0;    EY[e*Nverts+1] = y0+dy; EZ[e*Nverts+1] = z0+dz;
        EX[e*Nverts+2] = x0;    EY[e*Nverts+2] = y0;    EZ[e*Nverts+2] = z0+dz;
        EX[e*Nverts+3] = x0+dx; EY[e*Nverts+3] = y0+dy; EZ[e*Nverts+3] = z0+dz;

        elementInfo[e] = 1; // domain
        e++;

        //tet 5 (0,5,1,7)
        EToV[e*Nverts+0] = i0 + j0*NnX + k0*NnX*NnY;
        EToV[e*Nverts+1] = i1 + j0*NnX + k1*NnX*NnY;
        EToV[e*Nverts+2] = i1 + j0*NnX + k0*NnX*NnY;
        EToV[e*Nverts+3] = i1 + j1*NnX + k1*NnX*NnY;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;    EZ[e*Nverts+0] = z0;
        EX[e*Nverts+1] = x0+dx; EY[e*Nverts+1] = y0;    EZ[e*Nverts+1] = z0+dz;
        EX[e*Nverts+2] = x0+dx; EY[e*Nverts+2] = y0;    EZ[e*Nverts+2] = z0;
        EX[e*Nverts+3] = x0+dx; EY[e*Nverts+3] = y0+dy; EZ[e*Nverts+3] = z0+dz;

        elementInfo[e] = 1; // domain
        e++;

        //tet 6 (0,4,5,7)
        EToV[e*Nverts+0] = i0 + j0*NnX + k0*NnX*NnY;
        EToV[e*Nverts+1] = i0 + j0*NnX + k1*NnX*NnY;
        EToV[e*Nverts+2] = i1 + j0*NnX + k1*NnX*NnY;
        EToV[e*Nverts+3] = i1 + j1*NnX + k1*NnX*NnY;

        EX[e*Nverts+0] = x0;    EY[e*Nverts+0] = y0;    EZ[e*Nverts+0] = z0;
        EX[e*Nverts+1] = x0;    EY[e*Nverts+1] = y0;    EZ[e*Nverts+1] = z0+dz;
        EX[e*Nverts+2] = x0+dx; EY[e*Nverts+2] = y0;    EZ[e*Nverts+2] = z0+dz;
        EX[e*Nverts+3] = x0+dx; EY[e*Nverts+3] = y0+dy; EZ[e*Nverts+3] = z0+dz;

        elementInfo[e] = 1; // domain
        e++;
      }
    }
  }



  if (boundaryFlag != -1) { //-1 reserved for periodic case
    NboundaryFaces = 4*NX*NY + 4*NX*NZ + 4*NY*NZ;
    boundaryInfo = (hlong*) calloc(NboundaryFaces*(NfaceVertices+1), sizeof(hlong));

    hlong bcnt = 0;

    //top and bottom
    for(hlong j=0;j<NY;++j){
      for(hlong i=0;i<NX;++i){
        hlong vid1 = i + j*NnX +  0*NnX*NnY;
        hlong vid2 = i + j*NnX + NZ*NnX*NnY;

        boundaryInfo[bcnt*4+0] = boundaryFlag;
        boundaryInfo[bcnt*4+1] = vid1 + 0;
        boundaryInfo[bcnt*4+2] = vid1 + 1;
        boundaryInfo[bcnt*4+3] = vid1 + 1 + NnX;
        bcnt++;
        boundaryInfo[bcnt*4+0] = boundaryFlag;
        boundaryInfo[bcnt*4+1] = vid1 + 0;
        boundaryInfo[bcnt*4+2] = vid1 + 1 + NnX;
        boundaryInfo[bcnt*4+3] = vid1 + 0 + NnX;
        bcnt++;

        boundaryInfo[bcnt*4+0] = boundaryFlag;
        boundaryInfo[bcnt*4+1] = vid2 + 0;
        boundaryInfo[bcnt*4+2] = vid2 + 1;
        boundaryInfo[bcnt*4+3] = vid2 + 1 + NnX;
        bcnt++;
        boundaryInfo[bcnt*4+0] = boundaryFlag;
        boundaryInfo[bcnt*4+1] = vid2 + 0;
        boundaryInfo[bcnt*4+2] = vid2 + 1 + NnX;
        boundaryInfo[bcnt*4+3] = vid2 + 0 + NnX;
        bcnt++;
      }
    }
    //front and back
    for(hlong k=0;k<NZ;++k){
      for(hlong i=0;i<NX;++i){
        hlong vid1 = i +  0*NnX + k*NnX*NnY;
        hlong vid2 = i + NY*NnX + k*NnX*NnY;

        boundaryInfo[bcnt*4+0] = boundaryFlag;
        boundaryInfo[bcnt*4+1] = vid1 + 0;
        boundaryInfo[bcnt*4+2] = vid1 + 1;
        boundaryInfo[bcnt*4+3] = vid1 + 1 + NnX*NnY;
        bcnt++;
        boundaryInfo[bcnt*4+0] = boundaryFlag;
        boundaryInfo[bcnt*4+1] = vid1 + 0;
        boundaryInfo[bcnt*4+2] = vid1 + 1 + NnX*NnY;
        boundaryInfo[bcnt*4+3] = vid1 + 0 + NnX*NnY;
        bcnt++;

        boundaryInfo[bcnt*4+0] = boundaryFlag;
        boundaryInfo[bcnt*4+1] = vid2 + 0;
        boundaryInfo[bcnt*4+2] = vid2 + 1;
        boundaryInfo[bcnt*4+3] = vid2 + 1 + NnX*NnY;
        bcnt++;
        boundaryInfo[bcnt*4+0] = boundaryFlag;
        boundaryInfo[bcnt*4+1] = vid2 + 0;
        boundaryInfo[bcnt*4+2] = vid2 + 1 + NnX*NnY;
        boundaryInfo[bcnt*4+3] = vid2 + 0 + NnX*NnY;
        bcnt++;
      }
    }
    //left and right
    for(hlong k=0;k<NZ;++k){
      for(hlong j=0;j<NY;++j){
        hlong vid1 =  0 + j*NnX + k*NnX*NnY;
        hlong vid2 = NX + j*NnX + k*NnX*NnY;

        boundaryInfo[bcnt*4+0] = boundaryFlag;
        boundaryInfo[bcnt*4+1] = vid1 + 0*NnX;
        boundaryInfo[bcnt*4+2] = vid1 + 0*NnX + NnX*NnY;
        boundaryInfo[bcnt*4+3] = vid1 + 1*NnX + NnX*NnY;
        bcnt++;
        boundaryInfo[bcnt*4+0] = boundaryFlag;
        boundaryInfo[bcnt*4+1] = vid1 + 0*NnX;
        boundaryInfo[bcnt*4+2] = vid1 + 1*NnX + NnX*NnY;
        boundaryInfo[bcnt*4+3] = vid1 + 1*NnX;
        bcnt++;

        boundaryInfo[bcnt*4+0] = boundaryFlag;
        boundaryInfo[bcnt*4+1] = vid2 + 0*NnX;
        boundaryInfo[bcnt*4+2] = vid2 + 0*NnX + NnX*NnY;
        boundaryInfo[bcnt*4+3] = vid2 + 1*NnX + NnX*NnY;
        bcnt++;
        boundaryInfo[bcnt*4+0] = boundaryFlag;
        boundaryInfo[bcnt*4+1] = vid2 + 0*NnX;
        boundaryInfo[bcnt*4+2] = vid2 + 1*NnX + NnX*NnY;
        boundaryInfo[bcnt*4+3] = vid2 + 1*NnX;
        bcnt++;
      }
    }

  } else {
    NboundaryFaces = 0;
    boundaryInfo = NULL; // no boundaries
  }
}
