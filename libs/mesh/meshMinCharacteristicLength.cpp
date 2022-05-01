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

dfloat mesh_t::MinCharacteristicLength(){

  dfloat hmin = std::numeric_limits<dfloat>::max();
  for(dlong e=0;e<Nelements;++e){
    dfloat h = ElementCharacteristicLength(e);

    hmin = std::min(hmin, h);
  }

  // MPI_Allreduce to get global minimum h
  comm.Allreduce(hmin, Comm::Min);
  return hmin;
}

dfloat mesh_t::ElementCharacteristicLengthTri2D(dlong e) {

  dfloat h = std::numeric_limits<dfloat>::max();
  for(int f=0;f<Nfaces;++f){
    dlong sid = Nsgeo*(Nfaces*e + f);
    dfloat sJ   = sgeo[sid + SJID];
    dfloat invJ = sgeo[sid + IJID];

    // sJ = L/2, J = A/2,   sJ/J = L/A = L/(0.5*h*L) = 2/h
    // h = 2/(sJ/J)
    dfloat hest = 2.0/(sJ*invJ);

    h = std::min(h, hest);
  }
  return h;
}

dfloat mesh_t::ElementCharacteristicLengthQuad2D(dlong e) {

  dfloat h = std::numeric_limits<dfloat>::max();

  //sum weighted Jacobians to integrate over the element
  dfloat J = 0.0;
  for (int n=0;n<Np;n++)
    J += vgeo[Nvgeo*Np*e + n + Np*JWID];

  for(int f=0;f<Nfaces;++f){
    //sum weighted surface Jacobians to integrate over face
    dfloat sJ = 0.0;
    for (int i=0;i<Nfp;i++)
      sJ += sgeo[Nsgeo*(Nfaces*Nfp*e + Nfp*f + i) + WSJID];

    // sJ = L, J = A,   sJ/J = L/A = L/(h*L) = 1/h
    // h = 1/(sJ/J)
    dfloat hest = J/sJ;

    h = std::min(h, hest);
  }
  return h;
}

dfloat mesh_t::ElementCharacteristicLengthTet3D(dlong e) {

  dfloat h = std::numeric_limits<dfloat>::max();
  for(int f=0;f<Nfaces;++f){
    dlong sid = Nsgeo*(Nfaces*e + f);
    dfloat sJ   = sgeo[sid + SJID];
    dfloat invJ = sgeo[sid + IJID];

    // sJ = A/2, J = 3*V/4,   sJ/J = 2*A/3*V = 2*A/3*(A*h/3) = 2/h
    // h = 2/(sJ/J)
    dfloat hest = 2.0/(sJ*invJ);

    h = std::min(h, hest);
  }
  return h;
}

dfloat mesh_t::ElementCharacteristicLengthHex3D(dlong e) {

  dfloat h = std::numeric_limits<dfloat>::max();

  //sum weighted Jacobians to integrate over the element
  dfloat J = 0.0;
  for (int n=0;n<Np;n++)
    J += vgeo[Nvgeo*Np*e + n + Np*JWID];

  for(int f=0;f<Nfaces;++f){
    //sum weighted surface Jacobians to integrate over face
    dfloat sJ = 0.0;
    for (int i=0;i<Nfp;i++)
      sJ += sgeo[Nsgeo*(Nfaces*Nfp*e + Nfp*f + i) + WSJID];

    // sJ = L, J = A,   sJ/J = L/A = L/(h*L) = 1/h
    // h = 1/(sJ/J)
    dfloat hest = J/sJ;

    h = std::min(h, hest);
  }
  return h;
}

} //namespace libp
