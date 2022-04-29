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

#include "ins.hpp"

// compute RHS = beta*RHS + alpha*grad P
void ins_t::Gradient(const dfloat alpha, deviceMemory<dfloat>& o_P,
                     const dfloat beta,  deviceMemory<dfloat>& o_RHS,
                     const dfloat T){

  pTraceHalo.ExchangeStart(o_P, 1);

  // Compute Volume Contribution
  gradientVolumeKernel(mesh.Nelements,
                      mesh.o_vgeo,
                      mesh.o_D,
                      alpha,
                      beta,
                      o_P,
                      o_RHS);

  pTraceHalo.ExchangeFinish(o_P, 1);

  // Compute Surface Conribution
  gradientSurfaceKernel(mesh.Nelements,
                       mesh.o_sgeo,
                       mesh.o_LIFT,
                       mesh.o_vmapM,
                       mesh.o_vmapP,
                       mesh.o_EToB,
                       T,
                       mesh.o_x,
                       mesh.o_y,
                       mesh.o_z,
                       nu,
                       alpha,
                       o_P,
                       o_RHS);
}
