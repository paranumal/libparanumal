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

void mesh_t::SetupBoxLine1D(){

  // find a factorization size = size_x*size_y such that
  //  size_x>=size_y and are 'close' to one another
  int size_x = size;

  //determine (x) rank coordinates for this processes
  int rank_x=rank;

  //get global size from settings
  dlong NX;
  settings.getSetting("BOX GLOBAL NX", NX);

  //number of local elements in each dimension
  dlong nx;
  settings.getSetting("BOX NX", nx);

  if (NX <= 0) { //if the user hasn't given global sizes
    //set global size by multiplying local size by grid dims
    NX = nx * size_x;
    settings.changeSetting("BOX GLOBAL NX", std::to_string(NX));
  } else {
    //WARNING setting global sizes on input overrides any local sizes provided
    nx = NX/size_x + ((rank_x < (NX % size_x)) ? 1 : 0);
  }

  int boundaryFlag;
  settings.getSetting("BOX BOUNDARY FLAG", boundaryFlag);

  const int periodicFlag = (boundaryFlag == -1) ? 1 : 0;

  //grid physical sizes
  dfloat DIMX;
  settings.getSetting("BOX DIMX", DIMX);

  //element sizes
  dfloat dx = DIMX/NX;

  dlong offset_x = rank_x*(NX/size_x) + std::min(rank_x, (NX % size_x));

  //bottom corner of physical domain
  dfloat X0 = -DIMX/2.0 + offset_x*dx;

  //global number of nodes in each dimension
  hlong NnX = periodicFlag ? NX : NX+1; //lose a node when periodic (repeated node)

  // build an nx box grid
  Nnodes = NnX; //global node count
  Nelements = nx; //local

  EToV.malloc(Nelements*Nverts);
  EX.malloc(Nelements*Nverts);

  elementInfo.malloc(Nelements);

  #pragma omp parallel for
  for(int i=0;i<nx;++i){
    
    const dlong e = i;
    
    const hlong i0 = i+offset_x;
    const hlong i1 = (i+1+offset_x)%NnX;

    EToV[e*Nverts+0] = i0;
    EToV[e*Nverts+1] = i1;

    dfloat x0 = X0 + dx*i;

    dfloat *ex = EX.ptr()+e*Nverts;
    
    ex[0] = x0;   
    ex[1] = x0+dx;
    elementInfo[e] = 1; // domain
  }
  
  
  if (boundaryFlag != -1) { //-1 reserved for periodic case
    NboundaryFaces = 2;
    boundaryInfo.malloc(NboundaryFaces*(NfaceVertices+1));

    hlong bcnt = 0;

    //left and right
    hlong vid1 =  0;
    hlong vid2 = NX;

    boundaryInfo[bcnt*2+0] = boundaryFlag;
    boundaryInfo[bcnt*2+1] = vid1;
    bcnt++;
    
    boundaryInfo[bcnt*2+0] = boundaryFlag;
    boundaryInfo[bcnt*2+1] = vid2;
    bcnt++;
  }
  else{
    NboundaryFaces = 0; // no boundaries
  }
}

} //namespace libp
