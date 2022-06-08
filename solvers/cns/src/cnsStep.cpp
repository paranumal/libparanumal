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

#include "cns.hpp"

dfloat cns_t::MaxWaveSpeed(deviceMemory<dfloat>& o_Q, const dfloat T){

  //Note: if this is on the critical path in the future, we should pre-allocate this
  deviceMemory<dfloat> o_maxSpeed = platform.malloc<dfloat>(mesh.Nelements);

  maxWaveSpeedKernel(mesh.Nelements,
                     mesh.o_vgeo,
                     mesh.o_sgeo,
                     mesh.o_vmapM,
                     mesh.o_EToB,
                     gamma,
                     mu,
                     T,
                     mesh.o_x,
                     mesh.o_y,
                     mesh.o_z,
                     o_Q,
                     o_maxSpeed);

  const dfloat vmax = platform.linAlg().max(mesh.Nelements, o_maxSpeed, mesh.comm);

  return vmax;
}

//evaluate ODE rhs = f(q,t)
void cns_t::rhsf(deviceMemory<dfloat>& o_Q, deviceMemory<dfloat>& o_RHS, const dfloat T){

  // extract q trace halo and start exchange
  fieldTraceHalo.ExchangeStart(o_Q, 1);

  // compute volume contributions to gradients
  gradVolumeKernel(mesh.Nelements,
                   mesh.o_vgeo,
                   mesh.o_D,
                   o_Q,
                   o_gradq);

  // complete trace halo exchange
  fieldTraceHalo.ExchangeFinish(o_Q, 1);

  // compute surface contributions to gradients
  gradSurfaceKernel(mesh.Nelements,
                    mesh.o_sgeo,
                    mesh.o_LIFT,
                    mesh.o_vmapM,
                    mesh.o_vmapP,
                    mesh.o_EToB,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    T,
                    mu,
                    gamma,
                    o_Q,
                    o_gradq);

  // extract viscousStresses trace halo and start exchange
  gradTraceHalo.ExchangeStart(o_gradq, 1);

  // compute volume contribution to cns RHS
  if (cubature) {
    cubatureVolumeKernel(mesh.Nelements,
                         mesh.o_vgeo,
                         mesh.o_cubvgeo,
                         mesh.o_cubD,
                         mesh.o_cubPDT,
                         mesh.o_cubInterp,
                         mesh.o_cubProject,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         T,
                         mu,
                         gamma,
                         o_Q,
                         o_gradq,
                         o_RHS);
  } else {
    volumeKernel(mesh.Nelements,
                 mesh.o_vgeo,
                 mesh.o_D,
                 mesh.o_x,
                 mesh.o_y,
                 mesh.o_z,
                 T,
                 mu,
                 gamma,
                 o_Q,
                 o_gradq,
                 o_RHS);
  }

  // complete trace halo exchange
  gradTraceHalo.ExchangeFinish(o_gradq, 1);

  if (cubature) {
      cubatureSurfaceKernel(mesh.Nelements,
                            mesh.o_vgeo,
                            mesh.o_cubsgeo,
                            mesh.o_vmapM,
                            mesh.o_vmapP,
                            mesh.o_EToB,
                            mesh.o_intInterp,
                            mesh.o_intLIFT,
                            mesh.o_intx,
                            mesh.o_inty,
                            mesh.o_intz,
                            T,
                            mu,
                            gamma,
                            o_Q,
                            o_gradq,
                            o_RHS);
    } else {
      surfaceKernel(mesh.Nelements,
                    mesh.o_sgeo,
                    mesh.o_LIFT,
                    mesh.o_vmapM,
                    mesh.o_vmapP,
                    mesh.o_EToB,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    T,
                    mu,
                    gamma,
                    o_Q,
                    o_gradq,
                    o_RHS);
    }
}
