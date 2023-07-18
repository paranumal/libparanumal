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

// mesh is the local partition
int mesh_t::PointInclusionTest(dlong element, dfloat xs, dfloat ys, dfloat zs){

  if(elementType==Mesh::TRIANGLES){
    dlong id1 = element*Nverts + 0;
    dlong id2 = element*Nverts + 1;
    dlong id3 = element*Nverts + 2;
    dfloat x1 = EX[id1], x2 = EX[id2], x3 = EX[id3];
    dfloat y1 = EY[id1], y2 = EY[id2], y3 = EY[id3];

    // Xs = (1-r-s)*X1 + r*X2 + s*X3
    // [ (X2-X1)  (X3-X1) ]*[r;s] = Xs-X1
    // [r;s] = (1/J)*[ (y3-y1)  -(x3-x1); -(y2-y1)  (x2-x1) ]*[Xs-X1]
    // J = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)

    dfloat J =    (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1);
    dfloat re = ( (y3-y1)*(xs-x1) - (x3-x1)*(ys-y1))/J;
    dfloat se = (-(y2-y1)*(xs-x1) + (x2-x1)*(ys-y1))/J;

    if(0<=re && 0<=se && re+se<=1){
      return 1;
    }else{
      return 0;
    }
  }

  return 0;
}

}
