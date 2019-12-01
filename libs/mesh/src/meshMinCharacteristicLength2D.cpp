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

dfloat mesh2D::ElementCharacteristicLength(dlong e) {

  dfloat h = 1.e9;
  for(int f=0;f<Nfaces;++f){
    dlong sid = Nsgeo*(Nfaces*e + f);
    dfloat sJ   = sgeo[sid + SJID];
    dfloat invJ = sgeo[sid + IJID];

    // sJ = L/2, J = A/2,   sJ/J = L/A = L/(0.5*h*L) = 2/h
    // h = 0.5/(sJ/J)
    dfloat hest = 0.5/(sJ*invJ);

    h = mymin(h, hest);
  }
  return h;
}


dfloat mesh2D::MinCharacteristicLength(){

  dfloat hmin = 1e9;
  for(dlong e=0;e<Nelements;++e){
    dfloat h = ElementCharacteristicLength(e);

    hmin = mymin(hmin, h);
  }

  // MPI_Allreduce to get global minimum dt
  dfloat ghmin = 0.0;
  MPI_Allreduce(&hmin, &ghmin, 1, MPI_DFLOAT, MPI_MIN, comm);

  return ghmin;
}
