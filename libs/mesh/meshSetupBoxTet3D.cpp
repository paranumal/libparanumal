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

void mesh_t::SetupBoxTet3D(){

  // find a factorization size = size_x*size_y*size_z such that
  //  size_x>=size_y>=size_z are all 'close' to one another
  int size_x, size_y, size_z;
  Factor3(size, size_x, size_y, size_z);

  //determine (x,y,z) rank coordinates for this processes
  int rank_x=-1, rank_y=-1, rank_z=-1;
  RankDecomp3(size_x, size_y, size_z,
              rank_x, rank_y, rank_z,
              rank);

  //get global size from settings
  dlong NX, NY, NZ;
  settings.getSetting("BOX GLOBAL NX", NX);
  settings.getSetting("BOX GLOBAL NY", NY);
  settings.getSetting("BOX GLOBAL NZ", NZ);

  //number of local elements in each dimension
  dlong nx, ny, nz;
  settings.getSetting("BOX NX", nx);
  settings.getSetting("BOX NY", ny);
  settings.getSetting("BOX NZ", nz);

  if (NX*NY*NZ <= 0) { //if the user hasn't given global sizes
    //set global size by multiplying local size by grid dims
    NX = nx * size_x;
    NY = ny * size_y;
    NZ = nz * size_z;
    settings.changeSetting("BOX GLOBAL NX", std::to_string(NX));
    settings.changeSetting("BOX GLOBAL NY", std::to_string(NY));
    settings.changeSetting("BOX GLOBAL NZ", std::to_string(NZ));
  } else {
    //WARNING setting global sizes on input overrides any local sizes provided
    nx = NX/size_x + ((rank_x < (NX % size_x)) ? 1 : 0);
    ny = NY/size_y + ((rank_y < (NY % size_y)) ? 1 : 0);
    nz = NZ/size_z + ((rank_z < (NZ % size_z)) ? 1 : 0);
  }

  int boundaryFlag;
  settings.getSetting("BOX BOUNDARY FLAG", boundaryFlag);

  const int periodicFlag = (boundaryFlag == -1) ? 1 : 0;

  //grid physical sizes
  dfloat DIMX, DIMY, DIMZ;
  settings.getSetting("BOX DIMX", DIMX);
  settings.getSetting("BOX DIMY", DIMY);
  settings.getSetting("BOX DIMZ", DIMZ);

  //element sizes
  dfloat dx = DIMX/NX;
  dfloat dy = DIMY/NY;
  dfloat dz = DIMZ/NZ;

  dlong offset_x = rank_x*(NX/size_x) + std::min(rank_x, (NX % size_x));
  dlong offset_y = rank_y*(NY/size_y) + std::min(rank_y, (NY % size_y));
  dlong offset_z = rank_z*(NZ/size_z) + std::min(rank_z, (NZ % size_z));

  //bottom corner of physical domain
  dfloat X0 = -DIMX/2.0 + offset_x*dx;
  dfloat Y0 = -DIMY/2.0 + offset_y*dy;
  dfloat Z0 = -DIMZ/2.0 + offset_z*dz;

  //global number of nodes in each dimension
  hlong NnX = periodicFlag ? NX : NX+1; //lose a node when periodic (repeated node)
  hlong NnY = periodicFlag ? NY : NY+1; //lose a node when periodic (repeated node)
  hlong NnZ = periodicFlag ? NZ : NZ+1; //lose a node when periodic (repeated node)

  // build an nx x ny x nz box grid
  Nnodes = NnX*NnY*NnZ; //global node count
  Nelements = 6*nx*ny*nz; //local element count (each cube divided into 6 tets)

  EToV.malloc(Nelements*Nverts);
  EX.malloc(Nelements*Nverts);
  EY.malloc(Nelements*Nverts);
  EZ.malloc(Nelements*Nverts);

  elementInfo.malloc(Nelements);

  #pragma omp parallel for collapse(3)
  for(int k=0;k<nz;++k){
    for(int j=0;j<ny;++j){
      for(int i=0;i<nx;++i){

        dlong e = 6*(i + j*nx + k*nx*ny);

        const hlong i0 = i+offset_x;
        const hlong i1 = (i+1+offset_x)%NnX;
        const hlong j0 = j+offset_y;
        const hlong j1 = (j+1+offset_y)%NnY;
        const hlong k0 = k+offset_z;
        const hlong k1 = (k+1+offset_z)%NnZ;

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
    boundaryInfo.malloc(NboundaryFaces*(NfaceVertices+1));

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
    NboundaryFaces = 0; // no boundaries
  }
}

} //namespace libp
