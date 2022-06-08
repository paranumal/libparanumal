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

// compute RHS = beta*RHS + alpha*L(U)
void ins_t::Diffusion(const dfloat alpha, deviceMemory<dfloat>& o_U,
                      const dfloat beta,  deviceMemory<dfloat>& o_RHS,
                      const dfloat T) {

  //IPDG
  velocityGradientKernel(mesh.Nelements,
                        mesh.o_vgeo,
                        mesh.o_D,
                        o_U,
                        o_GU);

  // dfloat4 storage -> 4 entries
  vTraceHalo.ExchangeStart(o_GU, 4);

  if(mesh.NinternalElements)
    diffusionKernel(mesh.NinternalElements,
                   mesh.o_internalElementIds,
                   mesh.o_vgeo,
                   mesh.o_sgeo,
                   mesh.o_D,
                   mesh.o_LIFT,
                   mesh.o_vmapM,
                   mesh.o_vmapP,
                   mesh.o_EToB,
                   vTau,
                   T,
                   mesh.o_x,
                   mesh.o_y,
                   mesh.o_z,
                   nu,
                   alpha,
                   beta,
                   o_GU,
                   o_RHS);

  vTraceHalo.ExchangeFinish(o_GU, 4);

  if(mesh.NhaloElements)
    diffusionKernel(mesh.NhaloElements,
                    mesh.o_haloElementIds,
                    mesh.o_vgeo,
                    mesh.o_sgeo,
                    mesh.o_D,
                    mesh.o_LIFT,
                    mesh.o_vmapM,
                    mesh.o_vmapP,
                    mesh.o_EToB,
                    vTau,
                    T,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    nu,
                    alpha,
                    beta,
                    o_GU,
                    o_RHS);
}
