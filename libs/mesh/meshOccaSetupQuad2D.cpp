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

void meshQuad2D::OccaSetup(){

  this->mesh2D::OccaSetup();

  o_D = platform.malloc(Nq*Nq*sizeof(dfloat), D);

  o_S    = o_D; //dummy
  o_MM   = o_D; //dummy
  o_sM   = o_D; //dummy
  o_LIFT = o_D; //dummy

  o_vgeo = platform.malloc((Nelements+totalHaloPairs)*Nvgeo*Np*sizeof(dfloat), vgeo);
  o_sgeo = platform.malloc(Nelements*Nfaces*Nfp*Nsgeo*sizeof(dfloat), sgeo);
  o_ggeo = platform.malloc(Nelements*Np*Nggeo*sizeof(dfloat), ggeo);

  float *ggeo32 = (float*) calloc(Nelements*Np*Nggeo,sizeof(float));
  for(int n=0;n<Nelements*Np*Nggeo;++n){    
    ggeo32[n] = ggeo[n];
  }
  o_ggeo32 = platform.malloc(Nelements*Np*Nggeo*sizeof(float), ggeo32);
  free(ggeo32);
}
