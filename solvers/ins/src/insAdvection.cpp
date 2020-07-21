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

#include "ins.hpp"

// compute RHS = beta*RHS + alpha*N(U)
void ins_t::Advection(const dfloat alpha, occa::memory& o_U,
                      const dfloat beta,  occa::memory& o_RHS,
                      const dfloat T) {

  vTraceHalo->ExchangeStart(o_U, 1, ogs_dfloat);

  if (cubature)
    advectionVolumeKernel(mesh.Nelements,
                         mesh.o_vgeo,
                         mesh.o_cubvgeo,
                         mesh.o_cubD,
                         mesh.o_cubPDT,
                         mesh.o_cubInterp,
                         mesh.o_cubProject,
                         alpha,
                         beta,
                         o_U,
                         o_RHS);
  else
    advectionVolumeKernel(mesh.Nelements,
                         mesh.o_vgeo,
                         mesh.o_D,
                         alpha,
                         beta,
                         o_U,
                         o_RHS);

  vTraceHalo->ExchangeFinish(o_U, 1, ogs_dfloat);

  if (cubature)
    advectionSurfaceKernel(mesh.Nelements,
                          mesh.o_vgeo,
                          mesh.o_cubsgeo,
                          mesh.o_intInterp,
                          mesh.o_intLIFT,
                          mesh.o_vmapM,
                          mesh.o_vmapP,
                          mesh.o_EToB,
                          T,
                          mesh.o_intx,
                          mesh.o_inty,
                          mesh.o_intz,
                          nu,
                          alpha,
                          o_U,
                          o_RHS);
  else
    advectionSurfaceKernel(mesh.Nelements,
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
                          o_U,
                          o_RHS);
}
